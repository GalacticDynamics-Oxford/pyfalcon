/// Python interface for the fast-multipole gravity solver "falcON" (W.Dehnen, ApJL, 536, 9, 2000)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include "forces.h"  // from falcON

// compatibility with Python 3
#if PY_MAJOR_VERSION >= 3
#define PyInt_Check PyLong_Check
#endif

/// a convenience function for accessing an element of a PyArrayObject with the given data type
template<typename DataType>
inline DataType& pyArrayElem(void* arr, npy_intp ind)
{
    return *static_cast<DataType*>(PyArray_GETPTR1(static_cast<PyArrayObject*>(arr), ind));
}

/// same as above, but for a 2d array
template<typename DataType>
inline DataType& pyArrayElem(void* arr, npy_intp ind1, npy_intp ind2)
{
    return *static_cast<DataType*>(PyArray_GETPTR2(static_cast<PyArrayObject*>(arr), ind1, ind2));
}

/// the only public function provided by this module
static const char* docstring =
    "Compute the gravitational potential and acceleration created by the given set of particles.\n"
    "Arguments:\n"
    "pos - a 2d array of shape Nx3 containing particle positions.\n"
    "mass - a 1d array of length N containing individual particle masses, "
    "or a single number if all masses are identical. Note that masses cannot be zero.\n"
    "eps - a 1d array of length N containing individual softening lengths, "
    "or a single number if all softening lengths are identical.\n"
    "active (optional) - a 1d bool array of length N specifying 'active' particles: "
    "gravity is sourced by all particles, but accelerations and potential will be computed "
    "at the locations of active particles only. This does not seem to save CPU time, though. "
    "If not provided, all particles are considered active.\n"
    "theta (optional) - tree opening angle; smaller values provide better accuracy but "
    "are more computationally expensive. Default value is 0.6 and is suitable in most cases.\n"
    "kernel (optional) - choice of the softening kernel: 0 means Plummer softening, "
    "1 (default) is a recommended kernel with a steeper density cutoff (~r^-7), "
    "2 and 3 are non-positive-definite kernels with an even faster decay at large r.\n\n"
    "Return:\n"
    "a tuple of two arrays - accelerations (Kx3) and potential (K) for K active particles only.";

static PyObject *gravity(PyObject* /*self*/, PyObject* args, PyObject* named_args)
{
    static const char* keywords[] = {"pos", "mass", "eps", "active", "theta", "kernel", NULL};
    PyObject* pos_obj = NULL;
    PyObject* mass_obj= NULL;
    PyObject* eps_obj = NULL;
    PyObject* active_obj = NULL;
    double theta = falcON::Default::theta;
    int kernel = falcON::Default::kernel;
    if(!PyArg_ParseTupleAndKeywords(args, named_args, "OOO|Odi", const_cast<char **>(keywords),
        &pos_obj, &mass_obj, &eps_obj, &active_obj, &theta, &kernel))
    {
        //PyErr_SetString(PyExc_TypeError, "gravity: must provide pos, mass, eps");
        return NULL;
    }

    // pos must be a 2d array
    if(!PyArray_Check(pos_obj) ||
        PyArray_NDIM  ((PyArrayObject*)pos_obj) != 2 ||
        PyArray_DIM   ((PyArrayObject*)pos_obj, 1) != 3 ||
        !(PyArray_TYPE((PyArrayObject*)pos_obj) == NPY_FLOAT ||
        PyArray_TYPE  ((PyArrayObject*)pos_obj) == NPY_DOUBLE) )
    {
        PyErr_SetString(PyExc_TypeError, "Argument 'pos' must be a 2d array of shape Nx3");
        return NULL;
    }
    int pos_dtype = PyArray_TYPE((PyArrayObject*)pos_obj);
    npy_intp nbody = PyArray_DIM((PyArrayObject*)pos_obj, 0), nactive = nbody;
    if(sizeof(unsigned int)==4 && nbody > 2147483647) {
        PyErr_SetString(PyExc_TypeError, "Number of bodies may not exceed 2^31");
        return NULL;
    }

    // mass may be a single number or an array of the same length as pos
    double mass = 0;
    int mass_dtype = -1;
    if(PyFloat_Check(mass_obj) || PyInt_Check(mass_obj) || PyLong_Check(mass_obj)) {
        mass = PyFloat_AsDouble(mass_obj);
        if(PyErr_Occurred())
            return NULL;
    } else {
        if(!PyArray_Check(mass_obj) ||
            PyArray_NDIM  ((PyArrayObject*)mass_obj) != 1 ||
            PyArray_DIM   ((PyArrayObject*)mass_obj, 0) != nbody ||
            !(PyArray_TYPE((PyArrayObject*)mass_obj) == NPY_FLOAT ||
            PyArray_TYPE  ((PyArrayObject*)mass_obj) == NPY_DOUBLE) )
        {
            PyErr_SetString(PyExc_TypeError, "Argument 'mass' must be either a single number "
                "or a 1d array of the same length as 'pos'");
            return NULL;
        }
        mass_dtype = PyArray_TYPE((PyArrayObject*)mass_obj);
    }

    // eps may be a single number or an array of the same length as pos
    double eps = 0;
    int eps_dtype = -1;
    if(PyFloat_Check(eps_obj) || PyInt_Check(eps_obj) || PyLong_Check(eps_obj)) {
        eps = PyFloat_AsDouble(eps_obj);
        if(PyErr_Occurred())
            return NULL;
    } else {
        if(!PyArray_Check(eps_obj) ||
            PyArray_NDIM  ((PyArrayObject*)eps_obj) != 1 ||
            PyArray_DIM   ((PyArrayObject*)eps_obj, 0) != nbody ||
            !(PyArray_TYPE((PyArrayObject*)eps_obj) == NPY_FLOAT ||
            PyArray_TYPE  ((PyArrayObject*)eps_obj) == NPY_DOUBLE) )
        {
            PyErr_SetString(PyExc_TypeError, "Argument 'eps' must be either a single number "
                "or a 1d array of the same length as 'pos'");
            return NULL;
        }
        eps_dtype = PyArray_TYPE((PyArrayObject*)eps_obj);
    }

    // list of active particles may be provided (if not, all particles are considered active)
    if(active_obj &&
        (!PyArray_Check(active_obj) ||
        PyArray_NDIM((PyArrayObject*)active_obj) != 1 ||
        PyArray_DIM ((PyArrayObject*)active_obj, 0) != nbody ||
        PyArray_TYPE((PyArrayObject*)active_obj) != NPY_BOOL) )
    {
        PyErr_SetString(PyExc_TypeError, "Argument 'active', if provided, must be a 1d array "
            "of the same length as 'pos' with dtype=bool");
        return NULL;
    }

    // output arrays will be allocated once the number of active particles is known
    PyObject *acc_obj = NULL, *pot_obj = NULL;

    try{
        // falcON snapshot class holding particle attributes
        unsigned int nbodies[falcON::bodytype::NUM]={0, 0, (unsigned int)nbody};
        falcON::bodies bodies(nbodies, falcON::fieldset(  // which particle properties are allocated:
            falcON::fieldset::gravity |                   // pos,vel,mass,potential,acceleration,flags;
            (eps_dtype!=-1 ? falcON::fieldset::e : 0))    // if needed, also individual softening lengths
        );

        // load input arrays into bodies
        npy_intp i=0;
        for(falcON::body b=bodies.begin_all_bodies(); b; ++b,++i) {
            if(pos_dtype == NPY_FLOAT) {
                b.pos()[0] = pyArrayElem<float>(pos_obj, i, 0);
                b.pos()[1] = pyArrayElem<float>(pos_obj, i, 1);
                b.pos()[2] = pyArrayElem<float>(pos_obj, i, 2);
            } else { // pos_dtype == NPY_DOUBLE
                b.pos()[0] = pyArrayElem<double>(pos_obj, i, 0);
                b.pos()[1] = pyArrayElem<double>(pos_obj, i, 1);
                b.pos()[2] = pyArrayElem<double>(pos_obj, i, 2);
            }

            if(mass_dtype == -1)
                b.mass() = mass;
            else if(mass_dtype == NPY_FLOAT)
                b.mass() = pyArrayElem<float>(mass_obj, i);
            else // mass_dtype == NPY_DOUBLE
                b.mass() = pyArrayElem<double>(mass_obj, i);

            if(eps_dtype == NPY_FLOAT)
                b.eps() = pyArrayElem<float>(eps_obj, i);
            else if(eps_dtype == NPY_DOUBLE)
                b.eps() = pyArrayElem<double>(eps_obj, i);
            // otherwise do not use individual softening lengths

            if(active_obj == NULL || pyArrayElem<bool>(active_obj, i) )
                b.flag_as_active();
            else {
                b.unflag_active();
                nactive--;
            }
        }

        // allocate output arrays for active particles only
        npy_intp dims[2] = {nactive, 3};
        acc_obj = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        pot_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        if(!acc_obj || !pot_obj) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to allocate output arrays");
            Py_XDECREF(acc_obj);
            Py_XDECREF(pot_obj);
            return NULL;
        }

        // construct the octtree from all particles
        falcON::forces forces(&bodies,
            /*global softening length*/ eps,
            /*tree opening angle*/ theta,
            /*type of softening kernel*/ static_cast<falcON::kern_type>(kernel),
            /*whether to use individual softening lengths*/ eps_dtype!=-1);
        forces.grow();

        // compute the force for active particles only
        forces.approximate_gravity(/*use all particles if 'active' is not provided*/ active_obj==NULL);

        // store the accelerations and potential for active particles in output arrays
        i=0;
        for(falcON::body b=bodies.begin_active_bodies(); b; b.next_active(),++i) {
            pyArrayElem<double>(acc_obj, i, 0) = b.acc()[0];
            pyArrayElem<double>(acc_obj, i, 1) = b.acc()[1];
            pyArrayElem<double>(acc_obj, i, 2) = b.acc()[2];
            pyArrayElem<double>(pot_obj, i   ) = b.pot();
        }
        return Py_BuildValue("NN", acc_obj, pot_obj);
    }

    catch(std::exception& ex) {
        PyErr_SetString(PyExc_RuntimeError, ex.what());
        Py_XDECREF(acc_obj);
        Py_XDECREF(pot_obj);
        return NULL;
    }
}


static PyMethodDef methods[] = {
    { "gravity", (PyCFunction)gravity, METH_VARARGS | METH_KEYWORDS, docstring },
    { NULL, NULL, 0, NULL } };

// initialization code is different between Python 2 and Python 3
#if PY_MAJOR_VERSION < 3
// Python 2.6-2.7
typedef struct PyModuleDef {
    int m_base;
    const char* m_name;
    const char* m_doc;
    Py_ssize_t m_size;
    PyMethodDef *m_methods;
} PyModuleDef;
#define PyModuleDef_HEAD_INIT 0
#define PyModule_Create(def) Py_InitModule3((def)->m_name, (def)->m_methods, (def)->m_doc)
static PyObject* PyInit_pyfalcon(void);
PyMODINIT_FUNC initpyfalcon(void) { PyInit_pyfalcon(); }
static PyObject*
#else
// Python 3
PyMODINIT_FUNC
#endif
PyInit_pyfalcon() {
    static struct PyModuleDef moduledef = {PyModuleDef_HEAD_INIT, "pyfalcon",
         "Python interface for the falcON fast-multipole gravity solver", -1, methods};
    PyObject *module = PyModule_Create(&moduledef);
    import_array1(module);
    return module;
}
