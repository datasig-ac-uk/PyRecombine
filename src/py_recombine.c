//
// Created by sam on 16/10/24.
//
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_TARGET_VERSION NPY_2_0_API_VERSION
#define NPY_NO_DEPRECATED_API NPY_2_0_API_VERSION
#include "numpy/ndarrayobject.h"
// #include <numpy/npy_2_compat.h>

#if NPY_API_VERSION < NPY_2_0_API_VERSION
#error "This module requires NumPy 2.0 or later"
#endif

#include "_recombine.h"

PyDoc_STRVAR(py_recombine_doc,
             "recombine(ensemble, selector=(0,1,2,...no_points-1),"
             " weights = (1,1,..,1), degree = 1) ensemble is a numpy"
             " array of vectors of type NP_DOUBLE referred to as"
             " points, the selector is a list of indexes to rows in"
             " the ensemble, weights is a list of positive weights of"
             " equal length to the selector and defines an empirical"
             " measure on the points in the ensemble."
             " Returns (retained_indexes, new weights) The arrays"
             " index_array, weights_array are single index numpy arrays"
             " and must have the same dimension and represent the indexes"
             " of the vectors and a mass distribution of positive weights"
             " (and at least one must be strictly positive) on them."
             " The returned weights are strictly positive, have the"
             " same total mass - but are supported on a subset of the"
             " initial chosen set of locations. If degree has its default"
             " value of 1 then the vector data has the same integral under"
             " both weight distributions; if degree is k then both sets of"
             " weights will have the same moments for all degrees at most k;"
             " the indexes returned are a subset of indexes in input"
             " index_array and mass cannot be further recombined onto a"
             " proper subset while preserving the integral and moments."
             " The default is to index of all the points, the default"
             " weights are 1. on each point indexed."
             " The default degree is one."
);

static PyObject* py_recombine(PyObject* self, PyObject* args, PyObject* kwargs)
{
    // INTERNAL
    //
    int src_locations_built_internally = 0;
    int src_weights_built_internally = 0;
    // match the mean - or higher moments
    size_t stCubatureDegree;
    // max no points at end - computed below
    ptrdiff_t NoDimensionsToCubature;
    // parameter summaries
    ptrdiff_t no_locations;
    ptrdiff_t point_dimension;
    double total_mass = 0.;
    // the final output
    PyObject* out = NULL;

    // THE INPUTS
    // the data - a (large) enumeration of vectors obtained by making a list of vectors and converting it to an array.
    PyArrayObject* data;
    // a list of the rows of interest
    PyArrayObject* src_locations = NULL;
    // their associated weights
    PyArrayObject* src_weights = NULL;
    // match the mean - or higher moments
    ptrdiff_t CubatureDegree = 1;

    // Pre declaration of variables that are only used in the main branch.
    // The compiler is complaining when variables are declared and initialised
    // between the goto and label exit
    PyArrayObject* snk_locations = NULL;
    PyArrayObject* snk_weights = NULL;
    double* NewWeights;
    size_t* LOCATIONS;
    ptrdiff_t noKeptLocations;
    size_t* KeptLocations;

    // usage def recombine(array1, *args, degree=1)
    static const char* kwlist[] = {"ensemble", "selector", "weights", "degree", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!|O!O!n:recombine", (char**)kwlist, &PyArray_Type, &data,
                                     &PyArray_Type, &src_locations, &PyArray_Type, &src_weights, &CubatureDegree))
        return out;
    // DATA VALIDATION
    //
    if (data == NULL)
    {
        return NULL;
    }
    else if (PyArray_NDIM(data) != 2 || PyArray_DIM(data, 0) == 0 || PyArray_DIM(data, 1) == 0)
    {
        // present but badly formed
        PyErr_SetString(PyExc_ValueError, "data is badly formed");
        return NULL;
    }
    else if (src_locations != NULL && ((PyArray_NDIM(src_locations) != 1 || PyArray_DIM(src_locations, 0) == 0)))
    {
        // present but badly formed
        PyErr_SetString(PyExc_ValueError, "source locations badly formed");
        return NULL;
    }
    else if (src_weights != NULL && (PyArray_NDIM(src_weights) != 1 || PyArray_DIM(src_weights, 0) == 0))
    {
        // present but badly formed
        PyErr_SetString(PyExc_ValueError, "source weights badly formed");
        return NULL;
    }
    else if ((src_weights != NULL && src_locations != NULL) && !PyArray_SAMESHAPE(src_weights, src_locations))
    {
        // present well formed but of different length
        PyErr_SetString(PyExc_ValueError, "source weights and source locations have different shapes");
        return NULL;
    }
    else if (CubatureDegree < 1)
    {
        PyErr_SetString(PyExc_ValueError, "invalid cubature degree");
        return NULL;
    }
    stCubatureDegree = CubatureDegree; //(convert from signed to unsigned)
    // create default locations (ALL) if not specified
    if (src_locations == NULL)
    {
        npy_intp* d = PyArray_DIMS(data);
        //d[0] = PyArray_DIM(data, 0);
        src_locations = (PyArrayObject*)PyArray_SimpleNew(1, d, NPY_INTP);
        size_t* LOCS =  (size_t*) PyArray_DATA(src_locations);
        ptrdiff_t id;
        for (id = 0; id < d[0]; ++id)
            LOCS[id] = id;
        src_locations_built_internally = 1;
    }
    // create default weights (1. on each point) if not specified
    if (src_weights == NULL)
    {
        npy_intp d[1];
        d[0] = PyArray_DIM(src_locations, 0);
        src_weights = (PyArrayObject*)PyArray_SimpleNew(1, d, NPY_DOUBLE);
        double* WTS = (double*) PyArray_DATA(src_weights);
        ptrdiff_t id;
        for (id = 0; id < d[0]; ++id)
            WTS[id] = 1.;
        src_weights_built_internally = 1;
    }
    // make all data contiguous and type compliant (This only applies to external data - we know that our created arrays are fine
    // note this requires a deref at the end - so does the fill in of defaults - but we only do one or the other
    {
        data = (PyArrayObject*)PyArray_ContiguousFromObject((PyObject*)data, NPY_DOUBLE, 2, 2);
        if (!src_locations_built_internally)
            src_locations = (PyArrayObject*)PyArray_ContiguousFromObject((PyObject*)src_locations, NPY_INTP, 1, 1);
        if (!src_weights_built_internally)
            src_weights = (PyArrayObject*)PyArray_ContiguousFromObject((PyObject*)src_weights, NPY_DOUBLE, 1, 1);
    }


    // PREPARE INPUTS AS C ARRAYS
    ptrdiff_t no_datapoints = PyArray_DIM(data, 0);
    point_dimension = PyArray_DIM(data, 1);
    double* DATA = (double*) PyArray_DATA(data);


    LOCATIONS = (size_t*) PyArray_DATA(src_locations);
    double* WEIGHTS = (double*) PyArray_DATA(src_weights);

    // map locations from integer indexes to pointers to double
    no_locations = PyArray_DIM(src_locations, 0);
    double** LOCATIONS2 = (double**)malloc(no_locations * sizeof(double*));
    ptrdiff_t id;
    for (id = 0; id < no_locations; ++id)
    {
        // check for data out of range
        if (LOCATIONS[id] >= no_datapoints)
            goto exit;
        LOCATIONS2[id] = &DATA[LOCATIONS[id] * point_dimension];
    }
    // normalize the weights
    for (id = 0; id < no_locations; ++id)
        total_mass += WEIGHTS[id];
    for (id = 0; id < no_locations; ++id)
        WEIGHTS[id] /= total_mass;


    // NoDimensionsToCubature = the max number of points needed for cubature
    _recombineC(
        stCubatureDegree
        , point_dimension
        , 0 // tells _recombineC to return NoDimensionsToCubature the required buffer size
        , &NoDimensionsToCubature
        , NULL
        , NULL
        , NULL
        , NULL
    );
    // Prepare to call the reduction algorithm
    // a variable that will eventually be amended to to indicate the actual number of points returned
    noKeptLocations = NoDimensionsToCubature;
    // a buffer of size iNoDimensionsToCubature big enough to store array of indexes to the kept points
    KeptLocations = (size_t*)malloc(noKeptLocations * sizeof(size_t));
    // a buffer of size NoDimensionsToCubature to store the weights of the kept points
    NewWeights = (double*)malloc(noKeptLocations * sizeof(double));

    _recombineC(
        stCubatureDegree
        , point_dimension
        , no_locations
        , &noKeptLocations
        , (const void**)LOCATIONS2
        , WEIGHTS
        , KeptLocations
        , NewWeights
    );
    // un-normalise the weights
    for (id = 0; id < noKeptLocations; ++id)
        NewWeights[id] *= total_mass;
    // MOVE ANSWERS TO OUT
    // MAKE NEW OUTPUT OBJECTS
    npy_intp d[1];
    d[0] = noKeptLocations;

    snk_locations = (PyArrayObject*)PyArray_SimpleNew(1, d, NPY_INTP);
    snk_weights = (PyArrayObject*)PyArray_SimpleNew(1, d, NPY_DOUBLE);
    // MOVE OUTPUT FROM BUFFERS TO THESE OBJECTS
    memcpy(PyArray_DATA(snk_locations), KeptLocations, noKeptLocations * sizeof(size_t));
    memcpy(PyArray_DATA(snk_weights), NewWeights, noKeptLocations * sizeof(double));
    // RELEASE BUFFERS
    free(KeptLocations);
    free(NewWeights);
    // CREATE OUTPUT
    out = PyTuple_Pack(2, snk_locations, snk_weights);


exit:
    // CLEANUP
    free(LOCATIONS2);
    Py_DECREF(data);
    Py_DECREF(src_locations);
    Py_DECREF(src_weights);
    // EXIT
    return out;
    // USEFUL NUMPY EXAMPLES
    //https://stackoverflow.com/questions/56182259/how-does-one-acces-numpy-multidimensionnal-array-in-c-extensions/56233469#56233469
    //https://stackoverflow.com/questions/6994725/reversing-axis-in-numpy-array-using-c-api/6997311#6997311
    //https://stackoverflow.com/questions/6994725/reversing-axis-in-numpy-array-using-c-api/6997311#699731
}




static PyMethodDef py_recombine_methods[] = {
    { "recombine", (PyCFunction) py_recombine, METH_VARARGS | METH_KEYWORDS, py_recombine_doc},
    { NULL, NULL, 0, NULL }
};


static struct PyModuleDef py_recombine_module = {
    PyModuleDef_HEAD_INIT,
    "_recombine",
    "Recombine function for Python",
    -1,
    py_recombine_methods,
    NULL,
    NULL,
    NULL,
    NULL
};


PyMODINIT_FUNC PyInit__recombine(void) {

    PyObject* m;
    if (PyArray_ImportNumPyAPI() < 0) {
        return NULL;
    }

    m = PyModule_Create(&py_recombine_module);

    return m;
}