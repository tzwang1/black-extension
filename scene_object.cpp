/***********************************************************
        

        Starter code for Assignment 3
        

        Implements scene_object.h

***********************************************************/

#include "scene_object.h"
#include <cmath>

// Checks whether values for a, b, and c solve a quadratic equation
bool solveQuadratic(float &a, float &b, float &c, float &t){
    float discr = b*b - 4*a*c;
    float x0, x1;
    if(discr < 0){
        return false;
    } else if(discr == 0){
        x0 = -b / (2*a);
        x1 = -b / (2*a);
    } else {
        x0 = (-b + sqrt(discr))/(2*a);
        x1 = (-b - sqrt(discr))/(2*a);
    }

    t = fmin(x0, x1);
    if(t < 0){
        return false;
    }

    return true;
}

// Sets values for an intersecting point
Intersection set_intersect_values(float t, Vector3D normal, Point3D point, const Matrix4x4& worldToModel, const Matrix4x4& modelToWorld) {
    Intersection intersection;
    intersection.t_value = t;
    intersection.point = modelToWorld * point;
    normal = transNorm(worldToModel, normal);
    normal.normalize();
    intersection.normal = normal;
    intersection.none = false;
    return intersection;
}


bool UnitSquare::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

    // Transform the ray (origin, direction) to object space
    Point3D origin = worldToModel * ray.origin;
    Vector3D direction = worldToModel * ray.dir;

    // Since the normal of the unit square is (0, 0, 1), 
    // we only need to consider the z coordinate.
    double t =  -origin[2] / direction[2];

    // Invalid intersection
    if(t <= 0){
        return false;
    }

    // Calculate the point where the ray intersects the xy-plane.
    Point3D point = origin + t*direction;
    Vector3D normal = Vector3D(0,0,1);

    // Check if the point of intersection is within the unit square.
    if((-0.5 <= point[0]) && (point[0] <= 0.5) && (-0.5 <= point[1]) && (point[1] <= 0.5)){
        if (ray.intersection.none || t < ray.intersection.t_value) {
            // ray.intersection.point = modelToWorld * point;
            // normal = transNorm(worldToModel, normal);
            // normal.normalize();
            // ray.intersection.normal = normal;
            // ray.intersection.t_value = t;
            // ray.intersection.none = false;
            ray.intersection = set_intersect_values(t, normal, point, worldToModel, modelToWorld);
            return true;
        }
    }

    return false;
}


bool UnitSphere::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
    
    // coordinates for t if the ray intersects with the unit circle
    float t;
    // Transform the ray (origin, direction) to object space
    Point3D origin = worldToModel * ray.origin;
    Vector3D direction = worldToModel * ray.dir;

    Vector3D center = origin - Point3D(0,0,0);
    float a = direction.dot(direction);
    float b = 2 * direction.dot(center);
    float c = center.dot(center) - 1;

    if(!solveQuadratic(a, b, c, t)){
        return false;
    }

    // Intersection point
    Point3D point = origin + t*direction;

    Vector3D normal(point[0], point[1], point[2]);
    normal.normalize();

    if(ray.intersection.none || t < ray.intersection.t_value){
        // ray.intersection.t_value = t;
        // ray.intersection.point = modelToWorld * point;
        // normal = transNorm(worldToModel, normal);
        // normal.normalize();
        // ray.intersection.normal = normal;
        // ray.intersection.none = false;
        ray.intersection = set_intersect_values(t, normal, point, worldToModel, modelToWorld);
        return true;
    }
    
	return false;
}

bool intersects_top_or_bot(float t, Point3D origin, Vector3D direction, Intersection intersection){
    Point3D point = origin + t * direction;
    if(point[0] * point[0] + point[2] * point[2] < 1 && t > 0&& t < intersection.t_value){
        return true;
    }
    return false;
}


bool UnitCylinder::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld) {

    // coordinates for t if the ray intersects with the unit circle
    float t;
    
    // Transform the ray (origin, direction) to object space
    Point3D origin = worldToModel * ray.origin;
    Vector3D direction = worldToModel * ray.dir;

    float a = direction[0] * direction[0] + direction[2] * direction[2];
    float b = direction[0] * origin[0] + direction[2] * origin[2];
    float c = origin[0] * origin[0] + origin[2] * origin[2] - 1;

    if(!solveQuadratic(a, b, c, t)){
        return false;
    }

    // Intersection point
    Point3D point = origin + t*direction;

    Vector3D normal(point[0], 0, point[2]);
    normal.normalize();

    Intersection intersection;
    // Check if ray intersects with cylinder body
    if (point[1] < 0.5 && point[1] > -0.5){
        intersection = set_intersect_values(t, normal, point, worldToModel, modelToWorld);
    }

    // Check if ray intersects with the top or bottom part of the cylinder
    // float t_bot = (-0.5 - origin[1]) / direction[1];
    // Vector3D norm_bot = Vector3D(0, -1, 0);
    // if(intersects_top_or_bot(t_bot, origin, direction, intersection)){
    //     intersection = set_intersect_values(t_bot, normal, point, worldToModel, modelToWorld);
    // }

    // float t_top = (0.5 - origin[1]) / direction[1];
    // Vector3D norm_top = Vector3D(0, 1, 0);
    // if(intersects_top_or_bot(t_top, origin, direction, intersection)){
    //     intersection = set_intersect_values(t_top, normal, point, worldToModel, modelToWorld);
    // }

    if (!intersection.none && (ray.intersection.none || intersection.t_value < ray.intersection.t_value))
	{
		ray.intersection = intersection;
		return true;
	}

    return false;
}

void SceneNode::rotate(char axis, double angle) {
    Matrix4x4 rotation;
    double toRadian = 2 * M_PI / 360.0;
    int i;

    for (i = 0; i < 2; i++) {
        switch (axis) {
        case 'x':
            rotation[0][0] = 1;
            rotation[1][1] = cos(angle * toRadian);
            rotation[1][2] = -sin(angle * toRadian);
            rotation[2][1] = sin(angle * toRadian);
            rotation[2][2] = cos(angle * toRadian);
            rotation[3][3] = 1;
            break;
        case 'y':
            rotation[0][0] = cos(angle * toRadian);
            rotation[0][2] = sin(angle * toRadian);
            rotation[1][1] = 1;
            rotation[2][0] = -sin(angle * toRadian);
            rotation[2][2] = cos(angle * toRadian);
            rotation[3][3] = 1;
            break;
        case 'z':
            rotation[0][0] = cos(angle * toRadian);
            rotation[0][1] = -sin(angle * toRadian);
            rotation[1][0] = sin(angle * toRadian);
            rotation[1][1] = cos(angle * toRadian);
            rotation[2][2] = 1;
            rotation[3][3] = 1;
            break;
        }
        if (i == 0) {
            this->trans = this->trans * rotation;
            angle = -angle;
        } else {
            this->invtrans = rotation * this->invtrans;
        }
    }
}

void SceneNode::translate(Vector3D trans) {
    Matrix4x4 translation;

    translation[0][3] = trans[0];
    translation[1][3] = trans[1];
    translation[2][3] = trans[2];
    this->trans = this->trans * translation;
    translation[0][3] = -trans[0];
    translation[1][3] = -trans[1];
    translation[2][3] = -trans[2];
    this->invtrans = translation * this->invtrans;
}

void SceneNode::scale(Point3D origin, double factor[3]) {
    Matrix4x4 scale;

    scale[0][0] = factor[0];
    scale[0][3] = origin[0] - factor[0] * origin[0];
    scale[1][1] = factor[1];
    scale[1][3] = origin[1] - factor[1] * origin[1];
    scale[2][2] = factor[2];
    scale[2][3] = origin[2] - factor[2] * origin[2];
    this->trans = this->trans * scale;
    scale[0][0] = 1 / factor[0];
    scale[0][3] = origin[0] - 1 / factor[0] * origin[0];
    scale[1][1] = 1 / factor[1];
    scale[1][3] = origin[1] - 1 / factor[1] * origin[1];
    scale[2][2] = 1 / factor[2];
    scale[2][3] = origin[2] - 1 / factor[2] * origin[2];
    this->invtrans = scale * this->invtrans;
}
