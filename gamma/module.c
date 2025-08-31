#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include "gamma.h"


struct gpy_distribution {
    struct gamma_distribution dist; /* Base distribution object used by C */
    PyArrayObject            *data; /* Borrowed data reference */
};


struct gpy_results {
    struct gamma_results res;   /* The results buffer used by the C code */
    PyArrayObject       *arr;   /* NumPy array containing the gamma distrib. */
};


static bool gpy_get_long(PyObject *obj, const char *attr, long *res)
{
    PyObject *ptr;

    ptr = PyObject_GetAttrString(obj, attr);
    if (!ptr) {
        return false;
    }
    Py_DECREF(ptr);
    *res = PyLong_AsLong(ptr);
    return *res != -1 || !PyErr_Occurred();
}


static bool gpy_get_double(PyObject *obj, const char *attr, double *res)
{
    PyObject *ptr;

    ptr = PyObject_GetAttrString(obj, attr);
    if (!ptr) {
        return false;
    }
    Py_DECREF(ptr);
    *res = PyFloat_AsDouble(ptr);
    return *res != -1.0 || !PyErr_Occurred();
}


static bool gpy_get_bool(PyObject *obj, const char *attr, bool *res)
{
    PyObject *ptr;
    int value;

    ptr = PyObject_GetAttrString(obj, attr);
    if (!ptr) {
        return false;
    }
    Py_DECREF(ptr);
    value = PyObject_IsTrue(ptr);
    *res = value;
    return value != -1;
}


static bool gpy_get_arrayf(PyObject *obj, const char *attr,
                           size_t    len, double      arr[])
{
    PyObject *ptr;
    double *addr;
    size_t i;

    ptr = PyObject_GetAttrString(obj, attr);
    if (!ptr) {
        return false;
    }
    Py_DECREF(ptr);

    if (!PyArray_Check(ptr)) {
        PyErr_Format(PyExc_ValueError, "\"%s\" must be a NumPy array", attr);
        return false;
    }
    if (PyArray_TYPE((PyArrayObject *)ptr) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError,
                        "All NumPy arrays must have dtype=double");
        return false;
    }
    if ((size_t)PyArray_SIZE((PyArrayObject *)ptr) < len) {
        PyErr_Format(PyExc_IndexError,
                     "Expected %zu elements in NumPy array; found %zu",
                     len, PyArray_SIZE((PyArrayObject *)ptr));
        return false;
    }

    addr = PyArray_DATA((PyArrayObject *)ptr);
    for (i = 0; i < len; i++) {
        arr[i] = addr[i];
    }

    return true;
}


static bool gpy_load_norm(PyObject *obj, gamma_norm_t *norm)
{
    const char *value;
    PyObject *ptr;

    ptr = PyObject_GetAttrString(obj, "norm");
    if (!ptr) {
        return false;
    }
    Py_DECREF(ptr);

    value = PyUnicode_AsUTF8(ptr);
    if (!value) {
        return false;
    }

    if (!strcmp(value, "GLOBAL")) {
        *norm = GAMMA_NORM_GLOBAL;
    } else if (!strcmp(value, "LOCAL")) {
        *norm = GAMMA_NORM_LOCAL;
    } else if (!strcmp(value, "ABSOLUTE")) {
        *norm = GAMMA_NORM_ABSOLUTE;
    } else {
        PyErr_Format(PyExc_ValueError, "Normalization string \"%s\" is invalid",
                     value);
        return false;
    }
    return true;
}


static bool gpy_load_params(struct gamma_params *params, PyObject *obj)
{
    return gpy_get_double(obj, "diff", &params->diff)
        && gpy_get_double(obj, "dta", &params->dta)
        && gpy_get_double(obj, "threshold", &params->thrsh)
        && gpy_load_norm(obj, &params->norm)
        && gpy_get_bool(obj, "relative", &params->rel);
}


static bool gpy_load_options(struct gamma_options *opts, PyObject *obj)
{
    return gpy_get_bool(obj, "pass_only", &opts->pass_only)
        && gpy_get_long(obj, "pattern_shrinks", &opts->shrinks);
}


static bool gpy_load_distribution(struct gpy_distribution *dist, PyObject *obj)
{
    double buffer[9];
    npy_intp *dims;

    if (!gpy_get_arrayf(obj, "matrix", 9, buffer)) {
        return false;
    }
    dist->dist.matrix.cols[0].vec[0] = buffer[0];
    dist->dist.matrix.cols[0].vec[1] = buffer[1];
    dist->dist.matrix.cols[0].vec[2] = buffer[2];
    dist->dist.matrix.cols[0].vec[3] = 0.0;

    dist->dist.matrix.cols[1].vec[0] = buffer[3];
    dist->dist.matrix.cols[1].vec[1] = buffer[4];
    dist->dist.matrix.cols[1].vec[2] = buffer[5];
    dist->dist.matrix.cols[1].vec[3] = 0.0;

    dist->dist.matrix.cols[2].vec[0] = buffer[6];
    dist->dist.matrix.cols[2].vec[1] = buffer[7];
    dist->dist.matrix.cols[2].vec[2] = buffer[8];
    dist->dist.matrix.cols[2].vec[3] = 0.0;

    if (!gpy_get_arrayf(obj, "origin", 3, buffer)) {
        return false;
    }
    dist->dist.matrix.cols[3].vec[0] = buffer[0];
    dist->dist.matrix.cols[3].vec[1] = buffer[1];
    dist->dist.matrix.cols[3].vec[2] = buffer[2];
    dist->dist.matrix.cols[3].vec[3] = 1.0;

    if (!gpy_get_arrayf(obj, "spacing", 3, buffer)) {
        return false;
    }
    dist->dist.matrix.cols[0] = gamma_vec_muls(&dist->dist.matrix.cols[0], buffer[0]);
    dist->dist.matrix.cols[1] = gamma_vec_muls(&dist->dist.matrix.cols[1], buffer[1]);
    dist->dist.matrix.cols[2] = gamma_vec_muls(&dist->dist.matrix.cols[2], buffer[2]);

    dist->data = (PyArrayObject *)PyObject_GetAttrString(obj, "data");
    if (!dist->data) {
        return false;
    }
    Py_DECREF(dist->data);
    if (!PyArray_Check(dist->data)) {
        PyErr_SetString(PyExc_ValueError, "Dose data must be a NumPy array");
        return false;
    }

    if (PyArray_TYPE((PyArrayObject *)dist->data) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError,
                        "Dose data array must have dtype=double");
        return false;
    }
    dims = PyArray_DIMS((PyArrayObject *)dist->data);
    dist->dist.dims.idx[0] = dims[0];
    dist->dist.dims.idx[1] = dims[1];
    dist->dist.dims.idx[2] = dims[2];

    if (!gamma_distribution_set(&dist->dist, &dist->dist.matrix, &dist->dist.dims,
                                PyArray_DATA((PyArrayObject *)dist->data))) {
        PyErr_SetString(PyExc_ArithmeticError, "Affine matrix is singular");
        return false;
    }

    return true;
}


static bool gpy_write_long(long src, PyObject *obj, const char *attr)
{
    PyObject *node;
    bool res;

    node = PyLong_FromLong(src);
    if (!node) {
        return false;
    }
    res = !PyObject_SetAttrString(obj, attr, node);
    Py_DECREF(node);
    return res;
}


static bool gpy_write_double(double x, PyObject *obj, const char *attr)
{
    PyObject *node;
    bool res;

    node = PyFloat_FromDouble(x);
    if (!node) {
        return false;
    }
    res = !PyObject_SetAttrString(obj, attr, node);
    Py_DECREF(node);
    return res;
}


static bool gpy_write_array(PyArrayObject *arr, PyObject *obj, const char *attr)
{
    bool res;

    res = !PyObject_SetAttrString(obj, attr, (PyObject *)arr);
    Py_DECREF(arr);
    return res;
}


static bool gpy_write_results(struct gpy_results *res, PyObject *obj)
{
    return gpy_write_long(res->res.stats.total, obj, "total")
        && gpy_write_long(res->res.pass, obj, "passed")
        && gpy_write_double(res->res.stats.min, obj, "min")
        && gpy_write_double(res->res.stats.max, obj, "max")
        && gpy_write_double(res->res.stats.mean, obj, "mean")
        && gpy_write_double(res->res.stats.msqr, obj, "msqr")
        && gpy_write_array(res->arr, obj, "dist");
}


static PyObject *gpy_compute(PyObject *self, PyObject *args)
{
    PyObject *pyparms, *pyopts, *pyref, *pymeas, *pyres;
    struct gamma_params params;
    struct gamma_options opts;
    struct gpy_distribution ref, meas;
    struct gpy_results res;
    bool code;

    (void)self;
    if (!PyArg_ParseTuple(args, "OOOOO", &pyparms, &pyopts,
                                         &pyref, &pymeas, &pyres)) {
        return NULL;
    }

    if (!gpy_load_params(&params, pyparms)
     || !gpy_load_options(&opts, pyopts)
     || !gpy_load_distribution(&ref, pyref)
     || !gpy_load_distribution(&meas, pymeas)) {
        return NULL;
    }

    res.arr = (PyArrayObject *)PyArray_NewLikeArray(meas.data, NPY_KEEPORDER, NULL, 1);
    if (!res.arr) {
        return NULL;
    }
    res.res.dist = PyArray_DATA(res.arr);

    gamma_compute(&params, &opts, &ref.dist, &meas.dist, &res.res);
    code = gpy_write_results(&res, pyres);
    Py_DECREF(res.res.dist);

    return code ? Py_None : NULL;
}


PyMODINIT_FUNC PyInit_cgamma(void)
{
    static PyMethodDef methods[] = {
        {
            .ml_name  = "compute",
            .ml_meth  = gpy_compute,
            .ml_flags = METH_VARARGS,
            .ml_doc   = "Compute the gamma index",
        },
        { 0 }
    };
    static struct PyModuleDef module = {
        PyModuleDef_HEAD_INIT,
        .m_name    = "cgamma",
        .m_doc     = "Compute the gamma index",
        .m_methods = methods,
    };

    import_array();
    return PyModuleDef_Init(&module);
}
