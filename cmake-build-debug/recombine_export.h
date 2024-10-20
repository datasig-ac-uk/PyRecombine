
#ifndef RECOMBINE_EXPORT_H
#define RECOMBINE_EXPORT_H

#ifdef RECOMBINE_STATIC_DEFINE
#  define RECOMBINE_EXPORT
#  define RECOMBINE_NO_EXPORT
#else
#  ifndef RECOMBINE_EXPORT
#    ifdef _PyRecombine_EXPORTS
        /* We are building this library */
#      define RECOMBINE_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define RECOMBINE_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef RECOMBINE_NO_EXPORT
#    define RECOMBINE_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef RECOMBINE_DEPRECATED
#  define RECOMBINE_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef RECOMBINE_DEPRECATED_EXPORT
#  define RECOMBINE_DEPRECATED_EXPORT RECOMBINE_EXPORT RECOMBINE_DEPRECATED
#endif

#ifndef RECOMBINE_DEPRECATED_NO_EXPORT
#  define RECOMBINE_DEPRECATED_NO_EXPORT RECOMBINE_NO_EXPORT RECOMBINE_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef RECOMBINE_NO_DEPRECATED
#    define RECOMBINE_NO_DEPRECATED
#  endif
#endif

#endif /* RECOMBINE_EXPORT_H */
