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
    // the final output
    PyObject* out = NULL;

    PyObject* py_data = NULL;
    PyObject* py_locations = NULL;
    PyObject* py_weights = NULL;

    // THE INPUTS
    // the data - a (large) enumeration of vectors obtained by making a list of vectors and converting it to an array.
    PyArrayObject* data = NULL;
    // a list of the rows of interest
    PyArrayObject* src_locations = NULL;
    // their associated weights
    PyArrayObject* src_weights = NULL;
    // match the mean - or higher moments
    npy_intp CubatureDegree = 1;

    // Pre declaration of variables that are only used in the main branch.
    // The compiler is complaining when variables are declared and initialised
    // between the goto and label exit
    PyArrayObject* snk_locations = NULL;
    PyArrayObject* snk_weights = NULL;

    // usage def recombine(array1, *args, degree=1)
    static const char* kwlist[] = {"ensemble", "selector", "weights", "degree", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|OOn:recombine", (char**)kwlist,
            &py_data, &py_locations, &py_weights, &CubatureDegree)) {
        return out;
    }

    // DATA VALIDATION
    if (py_data == NULL) {
        PyErr_SetString(PyExc_ValueError, "data is required");
        return NULL;
    }

    // Get the easy check out of the way early
    if (CubatureDegree < 1) {
        PyErr_SetString(PyExc_ValueError, "invalid cubature degree");
        return NULL;
    }
    const size_t stCubatureDegree = CubatureDegree;


    // data = (PyArrayObject*) PyArray_FROM_OTF(py_data, NPY_DOUBLE,
    //     NPY_ARRAY_IN_ARRAY | NPY_ARRAY_ENSUREARRAY);
    // The input data must be a rank 2 array containing doubles
    data = (PyArrayObject*) PyArray_FromAny(
        py_data,
        PyArray_DescrFromType(NPY_DOUBLE), // desired type
        2, // min depth
        2, // max depth
        NPY_ARRAY_IN_ARRAY | NPY_ARRAY_ENSUREARRAY, // requirements
        NULL // context - unused
        );

    if (data == NULL) {
        // No references incremented if the data is not convertible to an array
        // and an error is set from the above function.
        return NULL;
    }

    if (PyArray_DIM(data, 0) == 0 || PyArray_DIM(data, 1) == 0) {
        PyErr_SetString(PyExc_ValueError, "data is empty");
        goto exit;
    }

    npy_intp no_datapoints = PyArray_DIM(data, 0);
    npy_intp point_dimension = PyArray_DIM(data, 1);
    npy_intp no_locations = no_datapoints;


    if (py_locations != NULL) {
        // Source locations provided, parse to array and check
        src_locations = (PyArrayObject*) PyArray_FromAny(
            py_locations,
            PyArray_DescrFromType(NPY_INTP), // desired type
            1,  // min depth
            1,  // max depth
            NPY_ARRAY_IN_ARRAY | NPY_ARRAY_ENSUREARRAY, // Requirements
            NULL // context - unused
            );

        if (src_locations == NULL) {
            // Error already set
            goto exit;
        }

        if (PyArray_DIM(src_locations, 0) == 0) {
            PyErr_SetString(PyExc_ValueError, "source locations is provided but empty");
            goto exit;
        }

        no_locations = PyArray_DIM(src_locations, 0);
        if (no_locations > no_datapoints) {
            PyErr_SetString(PyExc_ValueError, "source locations is provided but too large");
            goto exit;
        }

    } else {
        // If the src_locations array is empty then we should construct it ourselves.
        src_locations = (PyArrayObject*)PyArray_SimpleNew(1, &no_datapoints, NPY_INTP);
        npy_intp* LOCS =  (npy_intp*) PyArray_DATA(src_locations);

        for (npy_intp id = 0; id < no_datapoints; ++id) {
            LOCS[id] = id;
        }

        // src_locations_built_internally = 1;
    }


    double total_mass = 0.;
    if (py_weights != NULL) {
        // source weights provided, parse to array and check
        src_weights = (PyArrayObject*) PyArray_FromAny(
            py_weights,
            PyArray_DescrFromType(NPY_DOUBLE), // desired type
            1,  // min depth
            1,  // max depth
            NPY_ARRAY_ENSURECOPY | NPY_ARRAY_ENSUREARRAY, // Requirements
            NULL // context - unused
            );

        if (src_weights == NULL) {
            // Error already set
            goto exit;
        }

        if (PyArray_DIM(src_weights, 0) != no_locations) {
            PyErr_SetString(PyExc_ValueError,
                "source weights must match the number of locations");
            goto exit;
        }

        double* WTS = (double*) PyArray_DATA(src_weights);
        for (npy_intp id=0; id < no_locations; ++id) {
            total_mass += WTS[id];
        }

    } else {
        src_weights = (PyArrayObject*)PyArray_SimpleNew(1, &no_locations, NPY_DOUBLE);

        double* WTS = (double*) PyArray_DATA(src_weights);
        for (npy_intp id = 0; id < no_locations; ++id) {
            WTS[id] = 1.;
        }
        total_mass = (double) no_locations;

        // src_weights_built_internally = 1;
    }

    // PREPARE INPUTS AS C ARRAYS

    double* DATA = (double*) PyArray_DATA(data);

    size_t *LOCATIONS = (size_t *) PyArray_DATA(src_locations);
    double* WEIGHTS = (double*) PyArray_DATA(src_weights);

    // map locations from integer indexes to pointers to double
    double** LOCATIONS2 = (double**)malloc(no_locations * sizeof(double*));


    for (npy_intp id = 0; id < no_locations; ++id)
    {
        // check for data out of range
        if (LOCATIONS[id] >= no_locations) {
            PyErr_Format(PyExc_ValueError,
                "location %z out of range", LOCATIONS[id]);
            goto exit;
        }
        LOCATIONS2[id] = &DATA[LOCATIONS[id] * point_dimension];
    }

    // normalize the weights
    for (npy_intp id = 0; id < no_locations; ++id) {
        WEIGHTS[id] /= total_mass;
    }


    // NoDimensionsToCubature = the max number of points needed for cubature
    npy_intp NoDimensionsToCubature;
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
    // a variable that will eventually be amended to indicate the actual number of points returned
    npy_intp noKeptLocations = NoDimensionsToCubature;

    // a buffer of size iNoDimensionsToCubature big enough to store array of indexes to the kept points
    // size_t *KeptLocations = (size_t *) malloc(noKeptLocations * sizeof(size_t));

    // a buffer of size NoDimensionsToCubature to store the weights of the kept points
    // double* NewWeights = (double*)malloc(noKeptLocations * sizeof(double));

    snk_locations = (PyArrayObject*)PyArray_SimpleNew(1, &noKeptLocations, NPY_INTP);
    snk_weights = (PyArrayObject*)PyArray_SimpleNew(1, &noKeptLocations, NPY_DOUBLE);


    _recombineC(
        stCubatureDegree
        , point_dimension
        , no_locations
        , &noKeptLocations
        , (const void**)LOCATIONS2
        , WEIGHTS
        , (size_t*) PyArray_DATA(snk_locations)
        , (double*) PyArray_DATA(snk_weights)
    );



    // un-normalise the weights
    double* NewWeights = (double *) PyArray_DATA(snk_weights);
    for (npy_intp id = 0; id < noKeptLocations; ++id) {
        NewWeights[id] *= total_mass;
    }


    if (noKeptLocations != NoDimensionsToCubature) {
        PyArray_Dims dims = {&noKeptLocations, 1};
        if (PyArray_Resize(snk_locations, &dims, 0, NPY_CORDER) == NULL) {
            Py_XDECREF(snk_locations);
            Py_XDECREF(snk_weights);

            PyErr_SetString(PyExc_ValueError, "could not resize output arrays");
            goto exit;
        }
        if (PyArray_Resize(snk_weights, &dims, 0, NPY_CORDER) == NULL) {
            Py_XDECREF(snk_locations);
            Py_XDECREF(snk_weights);

            PyErr_SetString(PyExc_ValueError, "could not resize output arrays");
            goto exit;
        }

    }


    // memcpy(PyArray_DATA(snk_locations), KeptLocations, noKeptLocations * sizeof(size_t));
    // memcpy(PyArray_DATA(snk_weights), NewWeights, noKeptLocations * sizeof(double));

    // free(KeptLocations);
    // free(NewWeights);

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