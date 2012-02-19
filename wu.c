/* Based on code by Xiaolin Wu
 * Dept. of Computer Science
 * Univ. of Western Ontario
 * London, Ontario N6A 5B7
 * wu@csd.uwo.ca
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pam.h"

typedef struct color_cube {
    struct {int r,g,b;} min, max;
    double volume;
} color_cube;

typedef struct {
    float weight, r, g, b, moment;
} moments;

#define SIDE_SIZE 128

static moments m3d[SIDE_SIZE][SIDE_SIZE][SIDE_SIZE];

enum { DIRECTION_RED = 2,
       DIRECTION_GREEN = 1,
       DIRECTION_BLUE = 0,
};

// Converts the histogram to a series of moments.
static void calculate_moments()
{
    double area_w[SIDE_SIZE]={}, area_m[SIDE_SIZE]={},
           area_r[SIDE_SIZE]={}, area_g[SIDE_SIZE]={}, area_b[SIDE_SIZE]={};

    for (int red_idx = 1; red_idx < SIDE_SIZE; red_idx++) {
        for (int i = 0; i < SIDE_SIZE; i++) {
            area_w[i] = 0; area_r[i] = 0; area_g[i] = 0; area_b[i] = 0; area_m[i] = 0;
        }

        for (int green_idx = 1; green_idx < SIDE_SIZE; green_idx++) {
            double line_w = 0, line_r = 0, line_g = 0, line_b = 0, line_m = 0;

            for (int blue_idx = 1; blue_idx < SIDE_SIZE; blue_idx++) {
                area_w[blue_idx] += line_w += m3d[red_idx][green_idx][blue_idx].weight;
                area_r[blue_idx] += line_r += m3d[red_idx][green_idx][blue_idx].r;
                area_g[blue_idx] += line_g += m3d[red_idx][green_idx][blue_idx].g;
                area_b[blue_idx] += line_b += m3d[red_idx][green_idx][blue_idx].b;
                area_m[blue_idx] += line_m += m3d[red_idx][green_idx][blue_idx].moment;

                m3d[red_idx][green_idx][blue_idx].weight = m3d[red_idx-1][green_idx][blue_idx].weight + area_w[blue_idx];
                m3d[red_idx][green_idx][blue_idx].r      = m3d[red_idx-1][green_idx][blue_idx].r      + area_r[blue_idx];
                m3d[red_idx][green_idx][blue_idx].g      = m3d[red_idx-1][green_idx][blue_idx].g      + area_g[blue_idx];
                m3d[red_idx][green_idx][blue_idx].b      = m3d[red_idx-1][green_idx][blue_idx].b      + area_b[blue_idx];
                m3d[red_idx][green_idx][blue_idx].moment = m3d[red_idx-1][green_idx][blue_idx].moment + area_m[blue_idx];
            }
        }
    }
}

#define cube_volume_f(field) m3d[cube->max.r][cube->max.g][cube->max.b].field - \
                             m3d[cube->max.r][cube->max.g][cube->min.b].field - \
                             m3d[cube->max.r][cube->min.g][cube->max.b].field + \
                             m3d[cube->max.r][cube->min.g][cube->min.b].field - \
                             m3d[cube->min.r][cube->max.g][cube->max.b].field + \
                             m3d[cube->min.r][cube->max.g][cube->min.b].field + \
                             m3d[cube->min.r][cube->min.g][cube->max.b].field - \
                             m3d[cube->min.r][cube->min.g][cube->min.b].field

// Computes the volume of the cube in a specific moment.
static moments cube_volume(const color_cube *cube)
{
    return (moments) {
        .weight = cube_volume_f(weight),
        .r      = cube_volume_f(r),
        .g      = cube_volume_f(g),
        .b      = cube_volume_f(b),
        .moment = cube_volume_f(moment),
    };
}

#define cube_top_r(field) (m3d[position][cube->max.g][cube->max.b].field - \
                           m3d[position][cube->max.g][cube->min.b].field - \
                           m3d[position][cube->min.g][cube->max.b].field + \
                           m3d[position][cube->min.g][cube->min.b].field)

#define cube_top_g(field) (m3d[cube->max.r][position][cube->max.b].field - \
                           m3d[cube->max.r][position][cube->min.b].field - \
                           m3d[cube->min.r][position][cube->max.b].field + \
                           m3d[cube->min.r][position][cube->min.b].field)

#define cube_top_b(field) (m3d[cube->max.r][cube->max.g][position].field - \
                           m3d[cube->max.r][cube->min.g][position].field - \
                           m3d[cube->min.r][cube->max.g][position].field + \
                           m3d[cube->min.r][cube->min.g][position].field)

// Splits the cube in given position, and color direction.
static moments cube_top(const color_cube *cube, int direction, int position)
{
    switch (direction) {
        case DIRECTION_RED:
            return (moments) {
                .r = cube_top_r(r),
                .g = cube_top_r(g),
                .b = cube_top_r(b),
                .weight = cube_top_r(weight),
            };

        case DIRECTION_GREEN:
            return (moments) {
                .r = cube_top_g(r),
                .g = cube_top_g(g),
                .b = cube_top_g(b),
                .weight = cube_top_g(weight),
            };

        case DIRECTION_BLUE: default:
            return (moments) {
                .r = cube_top_b(r),
                .g = cube_top_b(g),
                .b = cube_top_b(b),
                .weight = cube_top_b(weight),
            };
    }
}

#define cube_bottom_r(field) (-m3d[cube->min.r][cube->max.g][cube->max.b].field + \
                               m3d[cube->min.r][cube->max.g][cube->min.b].field + \
                               m3d[cube->min.r][cube->min.g][cube->max.b].field - \
                               m3d[cube->min.r][cube->min.g][cube->min.b].field)

#define cube_bottom_g(field) (-m3d[cube->max.r][cube->min.g][cube->max.b].field + \
                               m3d[cube->max.r][cube->min.g][cube->min.b].field + \
                               m3d[cube->min.r][cube->min.g][cube->max.b].field - \
                               m3d[cube->min.r][cube->min.g][cube->min.b].field)

#define cube_bottom_b(field) (-m3d[cube->max.r][cube->max.g][cube->min.b].field + \
                               m3d[cube->max.r][cube->min.g][cube->min.b].field + \
                               m3d[cube->min.r][cube->max.g][cube->min.b].field - \
                               m3d[cube->min.r][cube->min.g][cube->min.b].field)

// Splits the cube in a given color direction at its minimum.
static moments cube_bottom(const color_cube *cube, int direction)
{
    switch (direction) {
        case DIRECTION_RED:
            return (moments) {
                .r = cube_bottom_r(r),
                .g = cube_bottom_r(g),
                .b = cube_bottom_r(b),
                .weight = cube_bottom_r(weight),
            };

        case DIRECTION_GREEN:
            return (moments) {
                .r = cube_bottom_g(r),
                .g = cube_bottom_g(g),
                .b = cube_bottom_g(b),
                .weight = cube_bottom_g(weight),
            };

        case DIRECTION_BLUE: default:
            return (moments) {
                .r = cube_bottom_b(r),
                .g = cube_bottom_b(g),
                .b = cube_bottom_b(b),
                .weight = cube_bottom_b(weight),
            };
    }
}

// Calculates statistical variance for a given cube.
static double variance(const color_cube *cube)
{
    if (cube->volume > 0) {
        const moments volume = cube_volume(cube);
        const double distance = volume.r*volume.r + volume.g*volume.g + volume.b*volume.b;
        return volume.moment - (distance/volume.weight);
    }
    return 0;
}

//	Finds the optimal (maximal) position for the cut.
static double maximal_position(const color_cube *cube, const moments *cube_vol, int direction, int first, int last, int *cut_p)
{
    moments bot = cube_bottom(cube, direction);

    double result = 0;
    *cut_p = -1;

    for (int position = first; position < last; position++) {
        moments top = cube_top(cube, direction, position);
        // determines the cube cut at a certain position
        double half_r = bot.r + top.r;
        double half_g = bot.g + top.g;
        double half_b = bot.b + top.b;
        double half_weight = bot.weight + top.weight;

        // the cube cannot be cut at bottom (this would lead to empty cube)
        if (half_weight > 0) {
            double half_dist = half_r*half_r + half_g*half_g + half_b*half_b;
            double temp = half_dist / half_weight;

            half_r = cube_vol->r - half_r;
            half_g = cube_vol->g - half_g;
            half_b = cube_vol->b - half_b;
            half_weight = cube_vol->weight - half_weight;

            if (half_weight > 0) {
                half_dist = half_r*half_r + half_g*half_g + half_b*half_b;
                temp += half_dist / half_weight;

                if (temp > result) {
                    result = temp;
                    *cut_p = position;
                }
            }
        }
    }

    return result;
}

// Cuts a cube with another one.
static int cut_cube(color_cube *first, color_cube *second)
{
    int direction, cut_r, cut_g, cut_b;

    moments cube_vol = cube_volume(first);

    double max_r = maximal_position(first, &cube_vol, DIRECTION_RED, first->min.r + 1, first->max.r, &cut_r);
    double max_g = maximal_position(first, &cube_vol, DIRECTION_GREEN, first->min.g + 1, first->max.g, &cut_g);
    double max_b = maximal_position(first, &cube_vol, DIRECTION_BLUE, first->min.b + 1, first->max.b, &cut_b);

    if (max_r >= max_g && max_r >= max_b) {
        direction = DIRECTION_RED;

        // cannot split empty cube
        if (cut_r < 0) return 0;
    }
    else if (max_g >= max_r && max_g >= max_b) {
        direction = DIRECTION_GREEN;
    }
    else {
        direction = DIRECTION_BLUE;
    }

    second->max.r = first->max.r;
    second->max.g = first->max.g;
    second->max.b = first->max.b;

    // cuts in a certain direction
    switch (direction) {
        case DIRECTION_RED:
            second->min.r = first->max.r = cut_r;
            second->min.g = first->min.g;
            second->min.b = first->min.b;
            break;

        case DIRECTION_GREEN:
            second->min.r = first->min.r;
            second->min.g = first->max.g = cut_g;
            second->min.b = first->min.b;
            break;

        case DIRECTION_BLUE: default:
            second->min.r = first->min.r;
            second->min.g = first->min.g;
            second->min.b = first->max.b = cut_b;
            break;
    }

    // determines the volumes after cut
    first->volume = (first->max.r - first->min.r)*(first->max.g - first->min.g)*(first->max.b - first->min.b);
    second->volume = (second->max.r - second->min.r)*(second->max.g - second->min.g)*(second->max.b - second->min.b);

    // the cut was successfull
    return 1;
}

// adds given color to histogram. multiplier is color's popularity.
static void add_color(rgb_pixel color, double multiplier)
{
    // element 0 is for base or marginal value
    int index_r = 1 + color.r * (SIDE_SIZE-1) / 256;
    int index_g = 1 + color.g * (SIDE_SIZE-1) / 256;
    int index_b = 1 + color.b * (SIDE_SIZE-1) / 256;

    m3d[index_r][index_g][index_b].weight += multiplier;
    m3d[index_r][index_g][index_b].r      += multiplier * color.r;
    m3d[index_r][index_g][index_b].g      += multiplier * color.g;
    m3d[index_r][index_g][index_b].b      += multiplier * color.b;
    m3d[index_r][index_g][index_b].moment += multiplier * (color.r*color.r + color.g*color.g + color.b*color.b);
}

static colormap* get_palette(int colors)
{
    color_cube cubes[colors];
    memset(cubes, 0, sizeof(cubes[0])*colors);

    // resets the reference maximums
    cubes[0].max.r = SIDE_SIZE-1;
    cubes[0].max.g = SIDE_SIZE-1;
    cubes[0].max.b = SIDE_SIZE-1;

    double volume_variance[colors];
    memset(volume_variance, 0, sizeof(volume_variance[0])*colors);

    int next = 0;

    // processes the cubes
    for (int cube_idx = 1; cube_idx < colors; cube_idx++) {
        // if cut is possible; make it
        if (cut_cube(&cubes[next], &cubes[cube_idx])) {
            volume_variance[next] = variance(&cubes[next]);
            volume_variance[cube_idx] = variance(&cubes[cube_idx]);
        }
        else { // the cut was not possible, revert the index
            volume_variance[next] = 0;
            cube_idx--;
        }

        next = 0;
        double max_variance = volume_variance[0];

        for (int i = 1; i <= cube_idx; i++) {
            if (volume_variance[i] > max_variance) {
                max_variance = volume_variance[i];
                next = i;
            }
        }

        if (max_variance <= 0.0) {
            colors = cube_idx + 1;
            break;
        }
    }


    colormap *result = pam_colormap(colors);

    // precalculates lookup tables
    for (int k = 0; k < colors; k++) {
        moments color = cube_volume(&cubes[k]);

        if (color.weight > 0) {
            rgb_pixel px = (rgb_pixel) {
                .r=(int)(color.r / color.weight),
                .g=(int)(color.g / color.weight),
                .b=(int)(color.b / color.weight),
                .a = 255,
            };
            result->palette[k].acolor = to_f(internal_gamma, px);
            result->palette[k].popularity = color.weight;
        }
    }

    return result;
}

colormap *wucut(histogram *hist, float min_opqaue_val, int colors)
{
    memset(m3d, 0, sizeof(m3d));

    for(int i=0; i < hist->size; i++) {
        if (hist->achv[i].acolor.a < min_opqaue_val) continue;
        add_color(to_rgb(internal_gamma, hist->achv[i].acolor), hist->achv[i].adjusted_weight);
    }

    // preprocess the colors
    calculate_moments();

    return get_palette(colors);
}
