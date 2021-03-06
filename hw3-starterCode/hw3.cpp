/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Yuka Murata
 * *************************
*/

#ifdef WIN32
    #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
    #include <GL/gl.h>
    #include <GL/glut.h>
#elif defined(__APPLE__)
    #include <OpenGL/gl.h>
    #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
    #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <iostream>

#include <vector>
#include <cmath>

using namespace std;

#define MAX_TRIANGLES 200000
#define MAX_SPHERES 1000
#define MAX_LIGHTS 1000

char * filename = NULL;

// different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480
double aspect_ratio = ((double) WIDTH) / ((double)HEIGHT);

// the field of view of the camera
#define fov 60.0

#define PI 3.14159265

// variables for soft shadow and antialiasing
#define NUM_ADDITIONAL_LIGHTS 16.0
#define NUM_SUBRAYS 2.0

// furthest solution t
#define MIN_T_MAX 1e8

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

struct Triangle
{
    Vertex v[3];
};

struct Sphere
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
};

struct Light
{
    double position[3];
    double color[3];
};

struct Color
{
	double r;
	double g;
	double b;

	Color () { r = g = b = 0; }

	Color (double r, double g, double b)
	{
		this->r = r;
		this->g = g;
		this->b = b;
	}

	Color add (Color color)
	{
		return Color(this->r + color.r, this->g + color.g, this->b + color.b);
	}

	Color mult (double scale)
	{
		return Color(this->r * scale, this->g * scale, this->b * scale);
	}

	Color mult (Color color)
	{
		return Color(this->r * color.r, this->g * color.g, this->b * color.b);
	}

	void clamp()
	{
		if (this->r > 1.0) this->r = 1.0;
		if (this->g > 1.0) this->g = 1.0;
		if (this->b > 1.0) this->b = 1.0;
	}
};

// Customized 3D vector
struct MyVector
{
	double x;
	double y;
	double z;

	// default constructor
    MyVector() { x = y = z = 0.0; }

	// customized constructor
    MyVector(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    // Add two vectors
    MyVector add(MyVector vec)
    {
        return MyVector(this->x + vec.x, this->y + vec.y, this->z + vec.z);
    }

    // Scale this vector
	MyVector mult(double scale)
    {
        return MyVector(scale * this->x, scale * this->y, scale * this->z);
    }

    // Negate this vector
    MyVector neg() { return this->mult(-1.0f); }

    // Subtract a vector from this vector
    MyVector minus(MyVector vec)
    {
    	return this->add(vec.neg());
    }

    // Get the length of this vector
    double magnitude(){
    	return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
    }

    // Normalize this vector
	MyVector normalize()
	{
		return this->mult(1.0f / this->magnitude());
	}

	// Dot product
	double dot(MyVector v)
	{
		return this->x * v.x + this->y * v.y + this->z * v.z;
	}

	// Cross product
	MyVector crossP(MyVector p)
	{
		double dx = this->y * p.z - this->z * p.y;
		double dy = - (this->x * p.z - this->z * p.x);
		double dz = this->x * p.y - this->y * p.x;

		return MyVector(dx, dy, dz);
	}
};

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// Variables
Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

// Get a direction vector from camera to the image plane
MyVector get_directionRay(double x, double y)
{
	MyVector direction;
	double degree = fov / 2 * PI / 180.0;

	// Set (x, y) = (0, 0) coordinate in the center of image plane
	double displacementX = - ((double) WIDTH) / ((double) HEIGHT) * tan(degree);
	double displacementY = - tan(degree);

	// calculate a direction ray
	double imagePlaneSize = 2 * tan(degree) / ((double)HEIGHT);
	direction.x = imagePlaneSize * ((double)x) + displacementX;
	direction.y = imagePlaneSize * ((double)y) + displacementY;
	direction.z = -1.0f;

	return direction.normalize();
}

// Get the minimum t value between two t values
double min(double t0, double t1)
{
	if (t0 < t1) return t0;
	else return t1;
}

MyVector get_vertex_position (int triangle_idx, int vertex_idx)
{
	return MyVector(triangles[triangle_idx].v[vertex_idx].position[0],
				  	triangles[triangle_idx].v[vertex_idx].position[1],
				  	triangles[triangle_idx].v[vertex_idx].position[2]);
}

void get_vertex_positions (int triangle_idx, MyVector &v0, MyVector &v1, MyVector &v2)
{
	v0 = get_vertex_position(triangle_idx, 0);
	v1 = get_vertex_position(triangle_idx, 1);
	v2 = get_vertex_position(triangle_idx, 2);
}

MyVector get_vertex_normal (int triangle_idx, int vertex_idx)
{
	return MyVector(triangles[triangle_idx].v[vertex_idx].normal[0],
				  	triangles[triangle_idx].v[vertex_idx].normal[1],
				  	triangles[triangle_idx].v[vertex_idx].normal[2]);
}

void get_vertex_normals (int triangle_idx, MyVector &n0, MyVector &n1, MyVector &n2)
{
	n0 = get_vertex_normal(triangle_idx, 0);
	n1 = get_vertex_normal(triangle_idx, 1);
	n2 = get_vertex_normal(triangle_idx, 2);
}

Color get_vertex_diffuse (int triangle_idx, int vertex_idx)
{
	return Color(triangles[triangle_idx].v[vertex_idx].color_diffuse[0],   // R
				 triangles[triangle_idx].v[vertex_idx].color_diffuse[1],   // G
				 triangles[triangle_idx].v[vertex_idx].color_diffuse[2]);  // B
}

void get_vertex_diffuses (int triangle_idx, Color &kd0, Color &kd1, Color &kd2)
{
	kd0 = get_vertex_diffuse(triangle_idx, 0);
	kd1 = get_vertex_diffuse(triangle_idx, 1);
	kd2 = get_vertex_diffuse(triangle_idx, 2);
}

Color get_vertex_specular (int triangle_idx, int vertex_idx)
{
	return Color(triangles[triangle_idx].v[vertex_idx].color_specular[0],   // R
				 triangles[triangle_idx].v[vertex_idx].color_specular[1],   // G
				 triangles[triangle_idx].v[vertex_idx].color_specular[2]);  // B
}

void get_vertex_speculars (int triangle_idx, Color &ks0, Color &ks1, Color &ks2)
{
	ks0 = get_vertex_specular(triangle_idx, 0);
	ks1 = get_vertex_specular(triangle_idx, 1);
	ks2 = get_vertex_specular(triangle_idx, 2);
}

Color get_sphere_diffuse (int sphere_idx)
{
	return Color(spheres[sphere_idx].color_diffuse[0],   // R
				 spheres[sphere_idx].color_diffuse[1],   // G
				 spheres[sphere_idx].color_diffuse[2]);  // B
}

void get_sphere_diffuses (int sphere_idx, Color &kd0, Color &kd1, Color &kd2)
{
	kd0 = get_sphere_diffuse(sphere_idx);
	kd1 = get_sphere_diffuse(sphere_idx);
	kd2 = get_sphere_diffuse(sphere_idx);
}

Color get_sphere_specular (int sphere_idx)
{
	return Color(spheres[sphere_idx].color_specular[0],   // R
				 spheres[sphere_idx].color_specular[1],   // G
				 spheres[sphere_idx].color_specular[2]);  // B
}

void get_sphere_speculars (int sphere_idx, Color &ks0, Color &ks1, Color &ks2)
{
	ks0 = get_sphere_specular(sphere_idx);
	ks1 = get_sphere_specular(sphere_idx);
	ks2 = get_sphere_specular(sphere_idx);
}

MyVector get_lightPosition (int light_idx)
{
	return MyVector(lights[light_idx].position[0],
					lights[light_idx].position[1],
					lights[light_idx].position[2]);
}

Color get_lightColor (int light_idx)
{
	return Color(lights[light_idx].color[0],
				 lights[light_idx].color[1],
				 lights[light_idx].color[2]);
}

// Get a shadow ray from the intersection point
MyVector get_shadowRay(MyVector intersect, int light_idx)
{
	// Shadow ray is (light position - intersection point)
	return get_lightPosition(light_idx).minus(intersect).normalize();
}

// Get normal of a triangle
MyVector get_triangle_normal(int triangle_idx)
{
	// Three verteces of a triangle
	MyVector p0, p1, p2;
	get_vertex_positions(triangle_idx, p0, p1, p2);

	MyVector p0p1 = p1.minus(p0);
	MyVector p0p2 = p2.minus(p0);
	
	return p0p1.crossP(p0p2);
}

// Get reflection vector
MyVector get_reflectionVector(MyVector l, MyVector n)
{
	// Evaluate r = 2(l.n)n - l
	return (n.mult(2 * l.dot(n)).minus(l)).normalize();
}

// Get interpolated value
Color interpolate (Color k0, Color k1, Color k2, double alpha, double beta, double gamma)
{
	// alpha * k0 + beta * k1 + gamma * k2
	return k0.mult(alpha).add(k1.mult(beta)).add(k2.mult(gamma));
}

// Inside test for a triangle
bool intersectTriangle(int triangle_idx, MyVector direction, double t, MyVector origin, double &alpha, double &beta)
{
	MyVector normal = get_triangle_normal(triangle_idx);
	double areaAll = normal.dot(normal);

	// Vertices of the triangle
	MyVector v0, v1, v2;
	get_vertex_positions(triangle_idx, v0, v1, v2);

	// Intersection point
	MyVector p = origin.add(direction.mult(t));

	// Inside test following

	// Edge between v0 and v1
    MyVector edge01 = v1.minus(v0);
    MyVector v0p = p.minus(v0);
    MyVector perpendicular = edge01.crossP(v0p);
    // out side of edge01
    if (normal.dot(perpendicular) < 0) return false;
 
 	// Edge between v1 and v2
    MyVector edge12 = v2.minus(v1);
    MyVector v1p = p.minus(v1);
    perpendicular = edge12.crossP(v1p);
    alpha = normal.dot(perpendicular);
    // out side of edge12
    if (alpha < 0) return false;
 
 	// Edge between v0 and v2
    MyVector edge20 = v0.minus(v2);
    MyVector v2p = p.minus(v2);
    perpendicular = edge20.crossP(v2p);
    beta = normal.dot(perpendicular);
    // out side of edge20
    if (beta < 0) return false;

    // Get alpha and beta from area ratios
    alpha /= areaAll;
    beta /= areaAll;
 
    return true;
}

// Inside test for a triangle
bool InsideTest(int triangle_idx, MyVector direction, double t, MyVector origin)
{
	// Placeholders
	double placeholder_alpha, placeholder_beta;

	return intersectTriangle(triangle_idx,
							 direction,
							 t,
							 origin,
							 placeholder_alpha,
							 placeholder_beta);
}

void get_sphere_info (int sphere_index, double &xc, double &yc, double &zc, double &radius)
{
	xc = spheres[sphere_index].position[0];
	yc = spheres[sphere_index].position[1];
    zc = spheres[sphere_index].position[2];
	radius = spheres[sphere_index].radius;
}

// Check if a ray is blocked by an object
bool blocked(MyVector shadowRay, MyVector intersection, int light_idx)
{
	// Threshold
	double epsilon = 1e-8;

	// Light position
	MyVector light = get_lightPosition(light_idx);

	// Distance between the light and intersection
	double dist_to_light = light.minus(intersection).magnitude();

	// Check intersection with any spheres
	for (int idx = 0; idx < num_spheres; idx++)
    {
    	// Sphere center position and radius
    	double xc, yc, zc, radius;
       	get_sphere_info(idx, xc, yc, zc, radius);

        double b = 2 * (shadowRay.x * (intersection.x - xc) +
        				shadowRay.y * (intersection.y - yc) +
        				shadowRay.z * (intersection.z - zc));
	    double c = pow(intersection.x - xc, 2.0) + pow(intersection.y - yc, 2.0) + pow(intersection.z - zc, 2.0) - pow(radius, 2.0);
       	double root = b * b - 4 * c;

	    // No solution for t
	   	if (root < 0) continue;

       	double t0 = (- b + sqrt(root)) / 2.0;
		double t1 = (- b - sqrt(root)) / 2.0;

		if (t0 > epsilon && shadowRay.mult(t0).magnitude() < dist_to_light) return true;
		if (t1 > epsilon && shadowRay.mult(t1).magnitude() < dist_to_light) return true;
    }

    // Check intersection with any triangles
    for (int triangle_idx = 0; triangle_idx < num_triangles; triangle_idx++)
    {
       	MyVector n = get_triangle_normal(triangle_idx).normalize();
       	double nd = n.dot(shadowRay);

       	// ray is parallel to the plane
       	if (-epsilon < nd && nd < epsilon) continue;

       	// Arbitrary point on the triangle in order to calculate d
        MyVector p0 = get_vertex_position(triangle_idx, 0);
        double d = n.dot(p0);

       	double t = (d - n.dot(intersection)) / nd;

       	// ray intersects a plane behind the scene
       	if (t <= epsilon) continue;
   		
   		// Inside test
       	if (t > epsilon &&
       		shadowRay.mult(t).magnitude() < dist_to_light &&
       		InsideTest(triangle_idx, shadowRay, t, intersection)) return true;
    }

    return false;
}

void add_randomLights(int  light_idx, int last_idx, double num_lights)
{
	double delta = 0.05;
	//int random = rand() % 5 + 1; // get random number between 1 - 5

	static int randomNumbers[] = { -4, -3, -2, -1, 1, 2, 3, 4} ;
	int index = rand() % (sizeof randomNumbers / sizeof *randomNumbers);
	int random = randomNumbers[index];

	lights[last_idx].position[0] = lights[light_idx].position[0] + random * delta;
	lights[last_idx].position[1] = lights[light_idx].position[1] + random * delta;
	lights[last_idx].position[2] = lights[light_idx].position[2] + random * delta;

	// evenly distribute light intensity
	lights[last_idx].color[0] = lights[light_idx].color[0] / num_lights;
	lights[last_idx].color[1] = lights[light_idx].color[1] / num_lights;
	lights[last_idx].color[2] = lights[light_idx].color[2] / num_lights;
}

// Generate lights around one light source for soft shadow
// Positions of additional lights are randomized (Much faster)
void addLights()
{
	int last_idx = num_lights;

	for (int light_idx = 0; light_idx < num_lights; light_idx++)
	{
		for (int i = 0; i < NUM_ADDITIONAL_LIGHTS; i++)
		{
			add_randomLights(light_idx, last_idx++, NUM_ADDITIONAL_LIGHTS);
		}

		// Original light source
		lights[light_idx].color[0] = lights[light_idx].color[0] / NUM_ADDITIONAL_LIGHTS;
		lights[light_idx].color[1] = lights[light_idx].color[1] / NUM_ADDITIONAL_LIGHTS;
		lights[light_idx].color[2] = lights[light_idx].color[2] / NUM_ADDITIONAL_LIGHTS;
	}

	num_lights = last_idx;
}

void check_sphere_intersection(MyVector direction, double &min_t_so_far, int &min_t_idx)
{
	double xc, yc, zc, radius, root, b, c;
	double t0, t1, min_t;

	for (int idx = 0; idx < num_spheres; idx++)
	{
		//  Sphere center position and radius
		get_sphere_info(idx, xc, yc, zc, radius);

		b = -2.0 * (direction.x * xc + direction.y * yc + direction.z * zc);
		c = pow(xc, 2.0) + pow(yc, 2.0) + pow(zc, 2.0) - pow(radius, 2.0);

		root = b * b - 4 * c;

		// When there is at least a solution for t
		if (root >= 0)
		{
			t0 = (- b + sqrt(root)) / 2.0;
			t1 = (- b - sqrt(root)) / 2.0;

			if (t0 > 0 && t1 > 0) min_t = min(t0, t1);
			else if (t0 > 0) min_t = t0;
			else if (t1 > 0) min_t = t1;
			else continue;

			if (min_t < min_t_so_far) 
			{
				// Save minimum t value and the sphere index at this time
				min_t_so_far = min_t;
				min_t_idx = idx;
			}
		}
	}
}

void check_triangle_intersection(MyVector direction, double &min_t_so_far, int &min_t_idx, bool &min_t_with_sphere, double &alpha, double &beta)
{
	double min_t;
	double this_alpha, this_beta;
	double epsilon = 1e-8;

	for (int triangle_idx = 0; triangle_idx < num_triangles; triangle_idx++)
	{
		MyVector n = get_triangle_normal(triangle_idx).normalize();
		double nd = n.dot(direction);

		if (nd < epsilon || epsilon < nd)
		{
			MyVector p0 = get_vertex_position(triangle_idx, 0);

			min_t = n.dot(p0) / nd;
			
			// ray intersects a plane including the triangle
			if (min_t > 0 && min_t < min_t_so_far)
			{
				if (intersectTriangle(triangle_idx, direction, min_t, MyVector(0, 0, 0), alpha, beta))
				{
					// Save minimum t value and the triangle index at the time
					// Also save alpha and beta values
					min_t_so_far = min_t;
					min_t_idx = triangle_idx;
					min_t_with_sphere = false;
					this_alpha = alpha;
					this_beta = beta;
				}
			}
		}
	}

	// Put alpha and beta values back
	alpha = this_alpha;
	beta = this_beta;
}

Color evaluate_sphere_phong(MyVector direction, MyVector shadowRay, int light_idx, double min_t, int min_t_idx)
{
	double xc, yc, zc, radius;
	MyVector normal, reflect;

	// update variables with final result
	get_sphere_info(min_t_idx, xc, yc, zc, radius);

	// Calculate surface normal
	normal = direction.mult(min_t).minus(MyVector(xc, yc, zc)).mult(1 / radius);

	// Get reflection
	reflect = get_reflectionVector(shadowRay, normal);

	//Evaluate local phong model
	double ln = shadowRay.dot(normal);
	if (ln < 0) ln = 0.0; // if l dot n is negative, make it 0

	Color diffuse = get_sphere_diffuse(min_t_idx).mult(ln);

	double rv = reflect.dot(direction.neg());
	if (rv < 0) rv = 0.0; // if r dot v is negative, make it 0
	Color specular = get_sphere_specular(min_t_idx).mult(pow(rv, spheres[min_t_idx].shininess));

	// Evaluate phong model
	Color current = get_lightColor(light_idx).mult(diffuse.add(specular));
	return current;//phong_light.add(current);
}

Color evaluate_triangle_phong(MyVector direction, MyVector shadowRay, int light_idx, int min_t_idx, double alpha, double beta)
{
	// Use interpolated normal
	MyVector normal0, normal1, normal2;
	get_vertex_normals(min_t_idx, normal0, normal1, normal2);

	// below means normal = alpha * normal0 + beta * normal1 + (1-alpha-beta) * normal2
	MyVector normal = normal0.mult(alpha).add(normal1.mult(beta)).add(normal2.mult(1-alpha-beta)).normalize();

	MyVector reflect = get_reflectionVector(shadowRay, normal);

	// kd's at each triangle vertex
	Color v0_kd, v1_kd, v2_kd;
	get_vertex_diffuses(min_t_idx, v0_kd, v1_kd, v2_kd);

	// Get diffuse component
	double ln = shadowRay.dot(normal);
	if (ln < 0) ln = 0; // if l dot n is negative, make it 0
	if (ln > 1) ln = 1.0;
	Color kd = interpolate(v0_kd, v1_kd, v2_kd, alpha, beta, 1-alpha-beta);
	Color diffuse = kd.mult(ln);;

	// ks' at each riangle vertex
	Color v0_ks, v1_ks, v2_ks;
	get_vertex_speculars(min_t_idx, v0_ks, v1_ks, v2_ks);

	// Get specular component
	double rv = reflect.dot(direction.neg());
	if (rv < 0) rv = 0.0; // if r dot v is negative, make it 0
	if (rv > 1) rv = 1.0;
	Color ks = interpolate(v0_ks, v1_ks, v2_ks, alpha, beta, 1-alpha-beta);

	double shininesses[3] = { triangles[min_t_idx].v[0].shininess,
								triangles[min_t_idx].v[1].shininess,
								triangles[min_t_idx].v[2].shininess };
	double interpolate_shi = alpha * shininesses[0] + beta * shininesses[1] + (1-alpha-beta) * shininesses[2];

	Color specular = ks.mult(pow(rv, interpolate_shi));

	// Evaluate phong model
	Color light = get_lightColor(light_idx);
	Color current = light.mult(diffuse.add(specular));
	return current;//phong_light.add(current);
}

/**
 * compute_phong_light
 *
 * @param direction - direction vector
 * @param min_t_with_sphere - closest intersection on sphere?
 * @param min_t - minimum t value
 * @param min_t_idx - sphere or triangle index when min_t
 * @param alpha - alpha value of barycentric coordinates
 * @param beta - beta value of barycentric coordinates
 *
 * Evaluate phong model at the closest intersection
 * by checking all the directions to light sources
 * Return Phong model for the ray
 */
Color compute_phong_light(MyVector direction, bool min_t_with_sphere, double min_t, int min_t_idx, double alpha, double beta)
{
	// local variables
	MyVector shadowRay;
	Color phong_light;

	// If a ray intersects with an object
	if (min_t != MIN_T_MAX)
	{
		// Check all possible light directions
		for (int light_idx = 0; light_idx < num_lights; light_idx++)
		{
			// shadow ray is a ray from intersection to a light
			shadowRay = get_shadowRay(direction.mult(min_t), light_idx);
	
			// When shadow ray is not blocked, evaluate the phong model
			if (!blocked(shadowRay, direction.mult(min_t), light_idx))
			{
				if (min_t_with_sphere) // Shading for a sphere
				{
					phong_light = phong_light.add(evaluate_sphere_phong(direction, shadowRay, light_idx, min_t, min_t_idx));
				}
				else // Shading for a triangle
				{
					phong_light = phong_light.add(evaluate_triangle_phong(direction, shadowRay, light_idx, min_t_idx, alpha, beta));
				}
			}
		}
	}
	else // Set white background
	{
		Color white = Color(1.0, 1.0, 1.0);
		phong_light = phong_light.add(white);
	}
	return phong_light;
}

/**
 * send_ray
 *
 * @param x - range of [0, WIDTH-1]
 * @param sub_x - range of [0, NUM_SUBRAYS-1]
 * @param y - range of [0, HEIGHT-1]
 * @param sub_y - range of [0, NUM_SUBRAYS-1]
 *
 * Shoot a ray and evaluate phong model
 * Return a phong light of the ray
 */
Color send_ray(int x, int sub_x, int y, int sub_y)
{
	// local variables
	double alpha, beta;
	double min_t;
	double min_t_so_far = MIN_T_MAX;
	int min_t_idx;
	bool min_t_with_sphere = true;

	// Step 1: Fire a ray from COP
	double del = 1 / NUM_SUBRAYS; // delta width between each subray
	MyVector direction = get_directionRay(x + sub_x * del, y + sub_y * del);

	// Step 2: Calculate closest intersection among objects
	check_sphere_intersection(direction, min_t_so_far, min_t_idx);
	check_triangle_intersection(direction, min_t_so_far, min_t_idx, min_t_with_sphere, alpha, beta);
	min_t = min_t_so_far; // update value

	// Step 3: For the closest intersection, evaluate phong model
	return compute_phong_light(direction, min_t_with_sphere, min_t, min_t_idx, alpha, beta);
}

/**
 * draw_scene
 * Shoot a ray to evaluate phong model and plot pixels
 */
void draw_scene()
{
	// Generate additional light sources for soft shadow
	addLights();

	// Iterate pixel by pixel
	for(unsigned int x=0; x<WIDTH; x++)
	{
		glPointSize(2.0);    
		glBegin(GL_POINTS);
		
		for(unsigned int y=0; y<HEIGHT; y++)
		{
			Color phong_lights; // Sum of sub phong lights

			// sub rays for antialiasing
			for (int sub_x = 0; sub_x < NUM_SUBRAYS; sub_x++)
			{
				for (int sub_y = 0; sub_y < NUM_SUBRAYS; sub_y++)
				{
					// Evaluate phong light of the ray
					Color phong_light = send_ray(x, sub_x, y, sub_y);

					// Sum up sub phong light to
					phong_lights = phong_lights.add(phong_light);
				}
			}

			// Fianl color = ave(diffuse + specular lights) + ambient light
			Color phong_sum = Color(phong_lights.r / pow(NUM_SUBRAYS, 2) + ambient_light[0],
									phong_lights.g / pow(NUM_SUBRAYS, 2) + ambient_light[1],
									phong_lights.b / pow(NUM_SUBRAYS, 2) + ambient_light[2]);
			phong_sum.clamp(); // clamp phong light in [0, 1]

			// Plot pixel with the phong light
			plot_pixel(x, y, phong_sum.r * 255, phong_sum.g * 255, phong_sum.b * 255);
		}

		glEnd();
		glFlush();
    }

    printf("Ray Tracying Completed!\n");
    fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
    glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    buffer[y][x][0] = r;
    buffer[y][x][1] = g;
    buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    plot_pixel_display(x,y,r,g,b);
    if(mode == MODE_JPEG)
        plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
    printf("Saving JPEG file: %s\n", filename);

    ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
    if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
        printf("Error in Saving\n");
    else 
        printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
    if(strcasecmp(expected,found))
    {
        printf("Expected '%s ' found '%s '\n", expected, found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check(check,str);
    fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
    printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check("rad:",str);
    fscanf(file,"%lf",r);
    printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
    char s[100];
    fscanf(file,"%s",s);
    parse_check("shi:",s);
    fscanf(file,"%lf",shi);
    printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
    FILE * file = fopen(argv,"r");
    int number_of_objects;
    char type[50];
    Triangle t;
    Sphere s;
    Light l;
    fscanf(file,"%i", &number_of_objects);

    printf("number of objects: %i\n",number_of_objects);

    parse_doubles(file,"amb:",ambient_light);

    for(int i=0; i<number_of_objects; i++)
    {
        fscanf(file,"%s\n",type);
        printf("%s\n",type);
        if(strcasecmp(type,"triangle")==0)
        {
            printf("found triangle\n");
            for(int j=0;j < 3;j++)
            {
                parse_doubles(file,"pos:",t.v[j].position);
                parse_doubles(file,"nor:",t.v[j].normal);
                parse_doubles(file,"dif:",t.v[j].color_diffuse);
                parse_doubles(file,"spe:",t.v[j].color_specular);
                parse_shi(file,&t.v[j].shininess);
            }

            if(num_triangles == MAX_TRIANGLES)
            {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            triangles[num_triangles++] = t;
        }
        else if(strcasecmp(type,"sphere")==0)
        {
            printf("found sphere\n");

            parse_doubles(file,"pos:",s.position);
            parse_rad(file,&s.radius);
            parse_doubles(file,"dif:",s.color_diffuse);
            parse_doubles(file,"spe:",s.color_specular);
            parse_shi(file,&s.shininess);

            if(num_spheres == MAX_SPHERES)
            {
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }
            spheres[num_spheres++] = s;
        }
        else if(strcasecmp(type,"light")==0)
        {
            printf("found light\n");
            parse_doubles(file,"pos:",l.position);
            parse_doubles(file,"col:",l.color);

            if(num_lights == MAX_LIGHTS)
            {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }
            lights[num_lights++] = l;
        }
        else
        {
            printf("unknown type in scene description:\n%s\n",type);
            exit(0);
        }
    }
    return 0;
}

void display() {}

void init()
{
    glMatrixMode(GL_PROJECTION);
    glOrtho(0,WIDTH,0,HEIGHT,1,-1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
    //hack to make it only draw once
    static int once=0;
    if(!once)
    {
        draw_scene();
        if(mode == MODE_JPEG)
            save_jpg();
    }
    once=1;
}

int main(int argc, char ** argv)
{
    if ((argc < 2) || (argc > 3))
    {    
        printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
        exit(0);
    }
    if(argc == 3)
    {
        mode = MODE_JPEG;
        filename = argv[2];
    }
    else if(argc == 2)
        mode = MODE_DISPLAY;

    glutInit(&argc,argv);
    loadScene(argv[1]);

    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(WIDTH,HEIGHT);
    int window = glutCreateWindow("Ray Tracer");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    init();
    glutMainLoop();
}
