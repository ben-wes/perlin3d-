#include <math.h>
#include <string.h>
#include "m_pd.h"

static t_class *perlin3d_tilde_class;

typedef struct _perlin3d_tilde {
    t_object  x_obj;
    t_outlet *x_out;
    int p[256];
} t_perlin3d_tilde;

unsigned int lcg_next(unsigned int *seed) {
    const unsigned int a = 1664525;
    const unsigned int c = 1013904223;
    *seed = (a * (*seed) + c); // Modulo 2^32 is implicit due to unsigned int overflow
    return *seed;
}

void shuffle(int *array, int n, unsigned int seed) {
    for (int i = n - 1; i > 0; i--) {
        unsigned int rand = lcg_next(&seed);
        int j = rand % (i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

void initPermutation(t_perlin3d_tilde *x, unsigned int seed) {
    int basePermutation[256];
    for (int i = 0; i < 256; i++) {
        basePermutation[i] = i;
    }
    shuffle(basePermutation, 256, seed);
    for (int i = 0; i < 256; i++) {
        x->p[i] = basePermutation[i];
    }
}

static inline int fastfloor(float x) {
    return x > 0 ? (int)x : (int)x - 1;
}

static inline uint8_t hash(t_perlin3d_tilde *x, int i) {
    return x->p[(uint8_t)i];
}

static inline t_float fade(t_float t) {
    return t * t * t * (t * (t * 6 - 15) + 10);
}

static inline t_float lerp(t_float t, t_float a, t_float b) {
    return a + t * (b - a);
}

static inline t_float grad(int hash, t_float x, t_float y, t_float z) {
    int h = hash & 15;
    t_float u = h < 8 ? x : y;
    t_float v = h < 4 ? y : h == 12 || h == 14 ? x : z;
    return ((h & 1) ? -u : u) + ((h & 2) == 0 ? v : -v);
}

t_float perlin3d(t_perlin3d_tilde *x, t_float xin, t_float yin, t_float zin) {
    // Using floorf for single-precision float, ensuring we stay within float precision
    int X = fastfloor(xin) & 255;
    int Y = fastfloor(yin) & 255;
    int Z = fastfloor(zin) & 255;

    xin -= fastfloor(xin);
    yin -= fastfloor(yin);
    zin -= fastfloor(zin);

    t_float u = fade(xin);
    t_float v = fade(yin);
    t_float w = fade(zin);

    int A = hash(x, X) + Y, AA = hash(x, A) + Z, AB = hash(x, A + 1) + Z;
    int B = hash(x, X + 1) + Y, BA = hash(x, B) + Z, BB = hash(x, B + 1) + Z;

    return lerp(w, lerp(v, lerp(u, grad(hash(x, AA), xin, yin, zin), grad(hash(x, BA), xin - 1, yin, zin)),
                        lerp(u, grad(hash(x, AB), xin, yin - 1, zin), grad(hash(x, BB), xin - 1, yin - 1, zin))),
                lerp(v, lerp(u, grad(hash(x, AA + 1), xin, yin, zin - 1), grad(hash(x, BA + 1), xin - 1, yin, zin - 1)),
                        lerp(u, grad(hash(x, AB + 1), xin, yin - 1, zin - 1), grad(hash(x, BB + 1), xin - 1, yin - 1, zin - 1))));
}

static t_int *perlin3d_tilde_perform(t_int *w) {
    t_perlin3d_tilde *x = (t_perlin3d_tilde *)(w[1]);
    t_float *in1 = (t_float *)(w[2]);
    t_float *in2 = (t_float *)(w[3]);
    t_float *in3 = (t_float *)(w[4]);
    t_sample *out = (t_sample *)(w[5]);
    int n = (int)(w[6]);

    while (n--) {
        t_float x_coord = *(in1++);
        t_float y_coord = *(in2++);
        t_float z_coord = *(in3++);
        *out++ = perlin3d(x, x_coord, y_coord, z_coord);
    }
    return (w + 7);
}

void perlin3d_tilde_dsp(t_perlin3d_tilde *x, t_signal **sp) {
    dsp_add(perlin3d_tilde_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}

void *perlin3d_tilde_new(void) {
    t_perlin3d_tilde *x = (t_perlin3d_tilde *)pd_new(perlin3d_tilde_class);
    initPermutation(x, 0);

    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);

    x->x_out = outlet_new(&x->x_obj, &s_signal);

    return (void *)x;
}

void perlin3d_tilde_setup(void) {
    perlin3d_tilde_class = class_new(gensym("perlin3d~"),
                                     (t_newmethod)perlin3d_tilde_new,
                                     0,
                                     sizeof(t_perlin3d_tilde),
                                     CLASS_DEFAULT,
                                     A_NULL);
    class_addmethod(perlin3d_tilde_class, (t_method)perlin3d_tilde_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(perlin3d_tilde_class, t_perlin3d_tilde, x_obj);
}
