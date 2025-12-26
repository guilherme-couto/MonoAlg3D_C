
/*
 * Octave3D.cpp
 *
 * Creates an 'octave noise' pattern (multiple octaves of basic Perlin
 * noise, evaluated at a set of points provided by the user. A grid with
 * unit spacing is used for the Perlin grid - input points should be
 * transformed appropriately to generate features of desired size and
 * orientation. The user may also specify the number of octaves, and the
 * factor successively applied to each octave. This seeded version of the
 * code also takes as input a permutation table of the numbers 0 to 255
 * and an m x 2 set of offsets in 2D space to apply to the Perlin grid for
 * each octave. m must therefore be at least as large as N_freq.
 *
 * Calling syntax in MATLAB is
 *
 *		noisefield = Octave2D(points, N_freq, fade_factor, Ps, offset)
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include "math.h"


//Define Perlin vectors (256 points distributed evenly around the unit sphere)
static double vecs_x[256] =
    {0.0513, -0.1160, -0.1404, 0.1475, 0.2168, -0.1424, -0.2865, 0.1119, 0.3471,
     -0.0614, -0.3955, -0.0050, 0.4289, 0.0832, -0.4450, -0.1690, 0.4423, 0.2583,
     -0.4200, -0.3468, 0.3782, 0.4302, -0.3178, -0.5045, 0.2403, 0.5660, -0.1482,
     -0.6116, 0.0444, 0.6385, 0.0676, -0.6449, -0.1837, 0.6296, 0.2999, -0.5921,
     -0.4118, 0.5329, 0.5152, -0.4534, -0.6061, 0.3556, 0.6808, -0.2422, -0.7361,
     0.1167, 0.7694, 0.0170, -0.7791, -0.1545, 0.7639, 0.2914, -0.7238, -0.4230,
     0.6593, 0.5447, -0.5721, -0.6523, 0.4645, 0.7419, -0.3395, -0.8101, 0.2008,
     0.8543, -0.0528, -0.8724, -0.1000, 0.8634, 0.2526, -0.8270, -0.4001, 0.7641,
     0.5377, -0.6761, -0.6608, 0.5654, 0.7652, -0.4354, -0.8472, 0.2899, 0.9039,
     -0.1332, -0.9332, -0.0296, 0.9338, 0.1935, -0.9053, -0.3534, 0.8484, 0.5040,
     -0.7644, -0.6405, 0.6559, 0.7585, -0.5259, -0.8540, 0.3785, 0.9238, -0.2182,
     -0.9655, 0.0499, 0.9777, 0.1212, -0.9596, -0.2896, 0.9118, 0.4501, -0.8354,
     -0.5975, 0.7329, 0.7272, -0.6072, -0.8349, 0.4622, 0.9173, -0.3025, -0.9717,
     0.1329, 0.9962, 0.0412, -0.9900, -0.2144, 0.9533, 0.3813, -0.8871, -0.5366,
     0.7935, 0.6755, -0.6755, -0.7935, 0.5366, 0.8871, -0.3813, -0.9533, 0.2144,
     0.9900, -0.0412, -0.9962, -0.1329, 0.9717, 0.3025, -0.9173, -0.4622, 0.8349,
     0.6072, -0.7272, -0.7329, 0.5975, 0.8354, -0.4501, -0.9118, 0.2896, 0.9596,
     -0.1212, -0.9777, -0.0499, 0.9655, 0.2182, -0.9238, -0.3785, 0.8540, 0.5259,
     -0.7585, -0.6559, 0.6405, 0.7644, -0.5040, -0.8484, 0.3534, 0.9053, -0.1935,
     -0.9338, 0.0296, 0.9332, 0.1332, -0.9039, -0.2899, 0.8472, 0.4354, -0.7652,
     -0.5654, 0.6608, 0.6761, -0.5377, -0.7641, 0.4001, 0.8270, -0.2526, -0.8634,
     0.1000, 0.8724, 0.0528, -0.8543, -0.2008, 0.8101, 0.3395, -0.7419, -0.4645,
     0.6523, 0.5721, -0.5447, -0.6593, 0.4230, 0.7238, -0.2914, -0.7639, 0.1545,
     0.7791, -0.0170, -0.7694, -0.1167, 0.7361, 0.2422, -0.6808, -0.3556, 0.6061,
     0.4534, -0.5152, -0.5329, 0.4118, 0.5921, -0.2999, -0.6296, 0.1837, 0.6449,
     -0.0676, -0.6385, -0.0444, 0.6116, 0.1482, -0.5660, -0.2403, 0.5045, 0.3178,
     -0.4302, -0.3782, 0.3468, 0.4200, -0.2583, -0.4423, 0.1690, 0.4450, -0.0832,
     -0.4289, 0.0050, 0.3955, 0.0614, -0.3471, -0.1119, 0.2865, 0.1424, -0.2168,
     -0.1475, 0.1404, 0.1160, -0.0513};

static double vecs_y[256] =
    {-0.0719, -0.0992, 0.1377, 0.1794, -0.1486, -0.2526, 0.1299, 0.3182, -0.0888,
     -0.3730, 0.0299, 0.4142, 0.0429, -0.4392, -0.1254, 0.4460, 0.2135, -0.4336,
     -0.3029, 0.4015, 0.3894, -0.3503, -0.4687, 0.2810, 0.5371, -0.1959, -0.5910,
     0.0975, 0.6275, 0.0108, -0.6444, -0.1254, 0.6400, 0.2421, -0.6136, -0.3567,
     0.5652, 0.4648, -0.4956, -0.5625, 0.4066, 0.6457, -0.3006, -0.7110, 0.1807,
     0.7556, -0.0506, -0.7773, -0.0855, 0.7746, 0.2233, -0.7469, -0.3581, 0.6945,
     0.4854, -0.6184, -0.6006, 0.5207, 0.6996, -0.4039, -0.7789, 0.2716, 0.8354,
     -0.1277, -0.8667, -0.0233, 0.8713, 0.1766, -0.8486, -0.3273, 0.7988, 0.4705,
     -0.7231, -0.6014, 0.6234, 0.7156, -0.5026, -0.8092, 0.3643, 0.8789, -0.2126,
     -0.9221, 0.0523, 0.9372, 0.1118, -0.9232, -0.2743, 0.8803, 0.4302, -0.8096,
     -0.5743, 0.7130, 0.7021, -0.5933, -0.8093, 0.4542, 0.8923, -0.2997, -0.9483,
     0.1347, 0.9754, 0.0357, -0.9724, -0.2061, 0.9394, 0.3712, -0.8770, -0.5257,
     0.7873, 0.6648, -0.6727, -0.7840, 0.5369, 0.8795, -0.3839, -0.9481, 0.2186,
     0.9877, -0.0461, -0.9969, -0.1283, 0.9754, 0.2990, -0.9238, -0.4607, 0.8436,
     0.6084, -0.7374, -0.7374, 0.6084, 0.8436, -0.4607, -0.9238, 0.2990, 0.9754,
     -0.1283, -0.9969, -0.0461, 0.9877, 0.2186, -0.9481, -0.3839, 0.8795, 0.5369,
     -0.7840, -0.6727, 0.6648, 0.7873, -0.5257, -0.8770, 0.3712, 0.9394, -0.2061,
     -0.9724, 0.0357, 0.9754, 0.1347, -0.9483, -0.2997, 0.8923, 0.4542, -0.8093,
     -0.5933, 0.7021, 0.7130, -0.5743, -0.8096, 0.4302, 0.8803, -0.2743, -0.9232,
     0.1118, 0.9372, 0.0523, -0.9221, -0.2126, 0.8789, 0.3643, -0.8092, -0.5026,
     0.7156, 0.6234, -0.6014, -0.7231, 0.4705, 0.7988, -0.3273, -0.8486, 0.1766,
     0.8713, -0.0233, -0.8667, -0.1277, 0.8354, 0.2716, -0.7789, -0.4039, 0.6996,
     0.5207, -0.6006, -0.6184, 0.4854, 0.6945, -0.3581, -0.7469, 0.2233, 0.7746,
     -0.0855, -0.7773, -0.0506, 0.7556, 0.1807, -0.7110, -0.3006, 0.6457, 0.4066,
     -0.5625, -0.4956, 0.4648, 0.5652, -0.3567, -0.6136, 0.2421, 0.6400, -0.1254,
     -0.6444, 0.0108, 0.6275, 0.0975, -0.5910, -0.1959, 0.5371, 0.2810, -0.4687,
     -0.3503, 0.3894, 0.4015, -0.3029, -0.4336, 0.2135, 0.4460, -0.1254, -0.4392,
     0.0429, 0.4142, 0.0299, -0.3730, -0.0888, 0.3182, 0.1299, -0.2526, -0.1486,
     0.1794, 0.1377, -0.0992, -0.0719};

static double vecs_z[256] =
    {-0.9961, -0.9883, -0.9805, -0.9727, -0.9648, -0.9570, -0.9492, -0.9414, -0.9336,
     -0.9258, -0.9180, -0.9102, -0.9023, -0.8945, -0.8867, -0.8789, -0.8711, -0.8633,
     -0.8555, -0.8477, -0.8398, -0.8320, -0.8242, -0.8164, -0.8086, -0.8008, -0.7930,
     -0.7852, -0.7773, -0.7695, -0.7617, -0.7539, -0.7461, -0.7383, -0.7305, -0.7227,
     -0.7148, -0.7070, -0.6992, -0.6914, -0.6836, -0.6758, -0.6680, -0.6602, -0.6523,
     -0.6445, -0.6367, -0.6289, -0.6211, -0.6133, -0.6055, -0.5977, -0.5898, -0.5820,
     -0.5742, -0.5664, -0.5586, -0.5508, -0.5430, -0.5352, -0.5273, -0.5195, -0.5117,
     -0.5039, -0.4961, -0.4883, -0.4805, -0.4727, -0.4648, -0.4570, -0.4492, -0.4414,
     -0.4336, -0.4258, -0.4180, -0.4102, -0.4023, -0.3945, -0.3867, -0.3789, -0.3711,
     -0.3633, -0.3555, -0.3477, -0.3398, -0.3320, -0.3242, -0.3164, -0.3086, -0.3008,
     -0.2930, -0.2852, -0.2773, -0.2695, -0.2617, -0.2539, -0.2461, -0.2383, -0.2305,
     -0.2227, -0.2148, -0.2070, -0.1992, -0.1914, -0.1836, -0.1758, -0.1680, -0.1602,
     -0.1523, -0.1445, -0.1367, -0.1289, -0.1211, -0.1133, -0.1055, -0.0977, -0.0898,
     -0.0820, -0.0742, -0.0664, -0.0586, -0.0508, -0.0430, -0.0352, -0.0273, -0.0195,
     -0.0117, -0.0039, 0.0039, 0.0117, 0.0195, 0.0273, 0.0352, 0.0430, 0.0508,
     0.0586, 0.0664, 0.0742, 0.0820, 0.0898, 0.0977, 0.1055, 0.1133, 0.1211,
     0.1289, 0.1367, 0.1445, 0.1523, 0.1602, 0.1680, 0.1758, 0.1836, 0.1914,
     0.1992, 0.2070, 0.2148, 0.2227, 0.2305, 0.2383, 0.2461, 0.2539, 0.2617,
     0.2695, 0.2773, 0.2852, 0.2930, 0.3008, 0.3086, 0.3164, 0.3242, 0.3320,
     0.3398, 0.3477, 0.3555, 0.3633, 0.3711, 0.3789, 0.3867, 0.3945, 0.4023,
     0.4102, 0.4180, 0.4258, 0.4336, 0.4414, 0.4492, 0.4570, 0.4648, 0.4727,
     0.4805, 0.4883, 0.4961, 0.5039, 0.5117, 0.5195, 0.5273, 0.5352, 0.5430,
     0.5508, 0.5586, 0.5664, 0.5742, 0.5820, 0.5898, 0.5977, 0.6055, 0.6133,
     0.6211, 0.6289, 0.6367, 0.6445, 0.6523, 0.6602, 0.6680, 0.6758, 0.6836,
     0.6914, 0.6992, 0.7070, 0.7148, 0.7227, 0.7305, 0.7383, 0.7461, 0.7539,
     0.7617, 0.7695, 0.7773, 0.7852, 0.7930, 0.8008, 0.8086, 0.8164, 0.8242,
     0.8320, 0.8398, 0.8477, 0.8555, 0.8633, 0.8711, 0.8789, 0.8867, 0.8945,
     0.9023, 0.9102, 0.9180, 0.9258, 0.9336, 0.9414, 0.9492, 0.9570, 0.9648,
     0.9727, 0.9805, 0.9883, 0.9961};


typedef struct Vector3_s {
   double x, y, z;
} Vector3;

/* Linear interpolation function. alpha is the position (alpha = 0 takes
   only starting component, alpha = 1 takes only ending component) */
inline double linearInterpolate(double fstart, double fend, double alpha)
{
    return fstart + alpha * (fend - fstart);
}

/* Smoothing function. Takes an input co-ordinate and smoothes it, in the
 * sense that proximity to 0 or 1 pushes the value closer to that extreme
 * limit. The function used is:
     x_smooth = 6 x^5 - 15 x^4 + 10 x^3                  */
inline double smoothInterpolate(double fstart, double fend, double alpha)
{
    double alpha_smoothed = alpha * alpha * alpha * (10 - alpha * (15 - 6 * alpha) );
    return linearInterpolate(fstart, fend, alpha_smoothed);
}

/* Dot product function. Takes the dot product of two 2D vectors */
inline double dotProduct3D(Vector3 v1, Vector3 v2)
{
    return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

/* Random vector selection function. Takes as input three integers, and uses
   hashing to convert these into a 'random' selection from 0-255 */
inline int selectVector(unsigned int x, unsigned int y, unsigned int z, int *P)
{

    // Hash with bitwise XOR
    return P[(x % 256)^P[(y % 256)^P[(z % 256)]]];

}

/* Noise calculation function. Takes an input point in 3D space and
 calculates the dot product with the four vectors surrounding this point */
double noise3D(double x, double y, double z, int *P)
{

    // Convert the input (x,y,z) into the corresponding integers
    int int_x = (int)floor(x);
    int int_y = (int)floor(y);
    int int_z = (int)floor(z);
    // Use these to find the "box co-ordinates (co-ords relative to bottom corner of box)
    double box_x = x - int_x;
    double box_y = y - int_y;
    double box_z = z - int_z;

    // Define vectors to the corners
    Vector3 v_dlo = {-box_x, -box_y, -box_z};
    Vector3 v_dli = {-box_x, -box_y, 1-box_z};
    Vector3 v_dro = {1-box_x, -box_y, -box_z};
    Vector3 v_dri = {1-box_x, -box_y, 1-box_z};
    Vector3 v_ulo = {-box_x, 1-box_y, -box_z};
    Vector3 v_uli = {-box_x, 1-box_y, 1-box_z};
    Vector3 v_uro = {1-box_x, 1-box_y, -box_z};
    Vector3 v_uri = {1-box_x, 1-box_y, 1-box_z};

    // Use hash function to assign gradient vectors to the corners
    unsigned int k_dlo = selectVector(int_x, int_y, int_z, P);
    unsigned int k_dli = selectVector(int_x, int_y, int_z+1, P);
    unsigned int k_dro = selectVector(int_x+1, int_y, int_z, P);
    unsigned int k_dri = selectVector(int_x+1, int_y, int_z+1, P);
    unsigned int k_ulo = selectVector(int_x, int_y+1, int_z, P);
    unsigned int k_uli = selectVector(int_x, int_y+1, int_z+1, P);
    unsigned int k_uro = selectVector(int_x+1, int_y+1, int_z, P);
    unsigned int k_uri = selectVector(int_x+1, int_y+1, int_z+1, P);

    Vector3 g_dlo = {vecs_x[k_dlo], vecs_y[k_dlo], vecs_z[k_dlo]};
    Vector3 g_dli = {vecs_x[k_dli], vecs_y[k_dli], vecs_z[k_dli]};
    Vector3 g_dro = {vecs_x[k_dro], vecs_y[k_dro], vecs_z[k_dro]};
    Vector3 g_dri = {vecs_x[k_dri], vecs_y[k_dri], vecs_z[k_dri]};
    Vector3 g_ulo = {vecs_x[k_ulo], vecs_y[k_ulo], vecs_z[k_ulo]};
    Vector3 g_uli = {vecs_x[k_uli], vecs_y[k_uli], vecs_z[k_uli]};
    Vector3 g_uro = {vecs_x[k_uro], vecs_y[k_uro], vecs_z[k_uro]};
    Vector3 g_uri = {vecs_x[k_uri], vecs_y[k_uri], vecs_z[k_uri]};

    // Calculate the values of the dot products
    double dlo_val = dotProduct3D(v_dlo, g_dlo);
    double dli_val = dotProduct3D(v_dli, g_dli);
    double dro_val = dotProduct3D(v_dro, g_dro);
    double dri_val = dotProduct3D(v_dri, g_dri);
    double ulo_val = dotProduct3D(v_ulo, g_ulo);
    double uli_val = dotProduct3D(v_uli, g_uli);
    double uro_val = dotProduct3D(v_uro, g_uro);
    double uri_val = dotProduct3D(v_uri, g_uri);

    // Perform linear interpolation of these dot products and return the result
    double bottom_outside_noise = smoothInterpolate(dlo_val, dro_val, box_x);
    double top_outside_noise = smoothInterpolate(ulo_val, uro_val, box_x);
    double bottom_inside_noise = smoothInterpolate(dli_val, dri_val, box_x);
    double top_inside_noise = smoothInterpolate(uli_val, uri_val, box_x);
    double outside_noise = smoothInterpolate(bottom_outside_noise, top_outside_noise, box_y);
    double inside_noise = smoothInterpolate(bottom_inside_noise, top_inside_noise, box_y);
    return smoothInterpolate(outside_noise, inside_noise, box_z);

}


/* Main runner function. Takes as input the points where Perlin noise is to
 * be evaluated, and calculates the value of octave noise at each point,
 * which is the combination of multiple 'octaves' of Perlin noise. The
 * result is scaled back to lie within [0,1] according to the theoretical
 * maximum and minimum values of 2D octave noise. */
void OctaveNoise3D(double *points, double *noisefield, mwSize n, unsigned int N_freq, double fade_factor, int *Ps, double *offsets)
{

    // Initialise loop variables
    mwSize i;
    unsigned int j;

    // Initialise multipliers
    int freq_mult;
    double scale_mult;

    // Initialise individual permutation tables
    int *P;

    // Calculate the scaling factor for normalisation
    // This is calculated using a geometric progression, with starting value sqrt(3)/2
    double normalisation_factor = 0.8660254 * (1 - (double)pow((double)fade_factor, (int)N_freq) ) / (1 - fade_factor);
    double normalisation_denominator = 0.5 / normalisation_factor;

    // Perform a loop over all input points
    for (i=0; i<n; i++) {

        // Initialise value of noisefield here as zero
        noisefield[i] = 0;

        // Loop over frequencies
        for (j=0; j<N_freq; j++) {

            // Set up multiplier for this frequency
            freq_mult = (int) pow((double)2 ,int(j));
            scale_mult = (double)pow((double)fade_factor, int(j));

            // Grab out this octave's permutation table (pointer pointing to shifted position along the full array Ps)
            P = &Ps[256*j];

            // Add value of noisefield at this location
            noisefield[i] = noisefield[i] + scale_mult*noise3D(points[3*i]*freq_mult - offsets[3*j], points[3*i+1]*freq_mult - offsets[3*j+1], points[3*i+2]*freq_mult - offsets[3*j+2], P);

        }

        // Convert this value of the noisefield to its [0,1] normalised equivalent
        noisefield[i] = ( noisefield[i] + normalisation_factor ) * normalisation_denominator;

    }
}

/* "Main" function used for interface between C and MATLAB */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    // Read in or otherwise initialise variables
    size_t n_points = mxGetN(prhs[0]);         // points supplied as rows
    double *points = mxGetPr(prhs[0]);
    unsigned int N_freq = mxGetScalar(prhs[1]);
    double fade_factor = mxGetScalar(prhs[2]);
    int *Ps = (int *)mxGetData(prhs[3]);
    double *offsets = mxGetPr(prhs[4]);
    double *noisefield;

    // Create the output matrix
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)n_points,mxREAL);

    // Create a pointer to this output matrix
    noisefield = mxGetPr(plhs[0]);

    // Call the PerlinNoise3D function, which will update the values of noisefield
    OctaveNoise3D(points, noisefield, (mwSize)n_points, N_freq, fade_factor, Ps, offsets);

}
