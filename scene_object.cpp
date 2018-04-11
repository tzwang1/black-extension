/***********************************************************
        
    Starter code for Assignment 3

    Implements scene_object.h

    Equations for primitives taken from
    https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html#SECTION00023200000000000000

***********************************************************/

#include "scene_object.h"
#include <cmath>

// Checks whether values for a, b, and c solve a quadratic equation
bool solveQuadratic(double &a, double &b, double &c, double &t) {
    double discr = b * b - 4 * a * c;
    double x0, x1;

    if (discr < 0) {
        return false;
    } else if (discr == 0) {
        x0 = -b / (2 * a);
        x1 = -b / (2 * a);
    } else {
        x0 = (-b + sqrt(discr)) / (2 * a);
        x1 = (-b - sqrt(discr)) / (2 * a);
    }
    // solve for parameter t by taking the minimum between the two solns
    // return false if the t is negative
    if (x0 < 0 && x1 < 0) {
        return false;
    } else if (x0 < 0) {
        t = x1;
    } else if (x1 < 0) {
        t = x0;
    } else {
        t = fmin(x0, x1);
    }

    return true;
}

// Sets values for an intersecting point
void set_intersect_values(Ray3D &ray, double t, Vector3D normal, Point3D point,
                          const Matrix4x4 &worldToModel, const Matrix4x4 &modelToWorld) {
    ray.intersection.t_value = t;
    ray.intersection.point = modelToWorld * point;
    normal = transNorm(worldToModel, normal);
    normal.normalize();
    ray.intersection.normal = normal;
    ray.intersection.none = false;
}

bool UnitSquare::intersect(Ray3D &ray, const Matrix4x4 &worldToModel,
                           const Matrix4x4 &modelToWorld) {
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
    double t = -origin[2] / direction[2];

    // Invalid intersection
    if (t <= 0) {
        return false;
    }

    // Calculate the point where the ray intersects the xy-plane.
    Point3D point = origin + t * direction;
    Vector3D normal(0.0, 0.0, 1.0);

    // Check if the point of intersection is within the unit square.
    if ((-0.5 <= point[0]) && (point[0] <= 0.5) && (-0.5 <= point[1]) && (point[1] <= 0.5)) {
        if (ray.intersection.none || t < ray.intersection.t_value) {
            set_intersect_values(ray, t, normal, point, worldToModel, modelToWorld);

            ray.intersection.tex_u = point[0] + 0.5;
            ray.intersection.tex_v = 0.5 - point[1];

            // std::cout << ray.intersection.tex_u << std::endl; 
            // std::cout << ray.intersection.tex_v<< std::endl;
            return true;
        }
    }

    return false;
}

bool UnitCube::intersect(Ray3D &ray, const Matrix4x4 &worldToModel,
                         const Matrix4x4 &modelToWorld) {
    // Intersect code for UnitCube, centered at the origin, with vertices at (1,0,0),
    // (0,1,0), (0,0,1), (0,0,0), (1,1,1), (1,1,0), (1,0,1), (0,1,1)
    // Transform the ray (origin, direction) to object spac
    Point3D origin = worldToModel * ray.origin;
    Vector3D direction = worldToModel * ray.dir;

    //Calculate all intersections with the cube
    double t_xmin = (0 - origin[0])/direction[0];
    double t_xmax = (1 - origin[0])/direction[0];

    double t_ymin = (0 - origin[1])/direction[1];
    double t_ymax = (1 - origin[1])/direction[1];

    double t_zmin = (0 - origin[2])/direction[2];
    double t_zmax = (1 - origin[2])/direction[2];

    //Swap the values if necessary
    if(t_xmin > t_xmax){
        std::swap(t_xmin, t_xmax);
    }

    if(t_ymin > t_ymax){
        std::swap(t_ymin, t_ymax);
    }
    
    if(t_zmin > t_zmax){
        std::swap(t_zmin, t_zmax);
    }

    // Return if the ray doesn't intersect
    if(t_xmin > t_ymax || t_ymin > t_xmax){
        return false;
    }
    double t_min;
    double t_max;
    if(t_ymin > t_xmin) {
        t_min = t_ymin;
    } else {
        t_min = t_xmin;
    }

    if(t_ymax < t_xmax) {
        t_max = t_ymax;
    } else {
        t_max = t_xmax;
    }

    // Return if the ray doesn't intersect
    if(t_min > t_zmax || t_zmin > t_max) {
        return false;
    }

    if(t_zmin > t_min) {
        t_min = t_zmin;
    }

    if(t_zmax < t_max) {
        t_max = t_zmax;
    }

    double t = t_min;
    if(t < 0){
        return false;
    }

    if(ray.intersection.none || t < ray.intersection.t_value){
        // std::cout << "Found intersection" << std::endl;
        Point3D point = origin + t * direction;
        Vector3D normal;
        if(abs(1 - int(point[0])) < 1e-3) {
            normal = Vector3D(1.0, 0.0, 0.0);
            ray.intersection.tex_u = point[1];
            ray.intersection.tex_v = point[2];
        // } else if(abs(0 - int(point[0])) < 1e-3) {
        //     normal = Vector3D(-1.0, 0.0, 0.0);
        //     // ray.intersection.tex_u = point[0];
        //     // ray.intersection.tex_v = point[0];
        } else if(abs(1 - int(point[1])) < 1e-3) {
            normal = Vector3D(0.0, 1.0, 0.0);
            ray.intersection.tex_u = point[0];
            ray.intersection.tex_v = point[2];
        // } else if(abs(0 - int(point[1])) < 1e-3) {
        //     normal = Vector3D(0.0, -1.0, 0.0);
        //     // ray.intersection.tex_u = point[0];
        //     // ray.intersection.tex_v = point[2];
        } else if(abs(1 - int(point[2])) < 1e-3) {
            normal = Vector3D(0.0, 0.0, 1.0);
            ray.intersection.tex_u = point[0];
            ray.intersection.tex_v = point[1];
        // } else if(abs(0 - int(point[2])) < 1e-3) {
        //     normal = Vector3D(0.0, 0.0, -1.0);
        //     // ray.intersection.tex_u = point[0];
        //     // ray.intersection.tex_v = point[1];
        } else {
            // std::cout << "Something went wrong... this is what the point looks like" << std::endl;
            // std::cout << point << std::endl;
            // return false;
        }
        set_intersect_values(ray, t, normal, point, worldToModel, modelToWorld);
        return true;
    }  
    return false;
}

bool UnitSphere::intersect(Ray3D &ray, const Matrix4x4 &worldToModel,
                           const Matrix4x4 &modelToWorld) {
    // TODO: implement intersection code for UnitSphere, which is centred
    // on the origin.
    //
    // Your goal here is to fill ray.intersection with correct values
    // should an intersection occur.  This includes intersection.point,
    // intersection.normal, intersection.none, intersection.t_value.
    //
    // HINT: Remember to first transform the ray into object space
    // to simplify the intersection test.

    // t is the parameter for the point that intersects the sphere
    double t;
    // Transform the ray (origin, direction) to object space
    Point3D origin = worldToModel * ray.origin;
    Vector3D direction = worldToModel * ray.dir;

    Vector3D center = origin - Point3D(0, 0, 0);
    double a = direction.dot(direction);
    double b = 2 * direction.dot(center);
    double c = center.dot(center) - 1;

    if (!solveQuadratic(a, b, c, t)) {
        return false;
    }

    // Intersection point
    Point3D point = origin + t * direction;

    Vector3D normal(point[0], point[1], point[2]);
    normal.normalize();

    if (ray.intersection.none || t < ray.intersection.t_value) {
        set_intersect_values(ray, t, normal, point, worldToModel, modelToWorld);
        Vector3D n = ray.intersection.normal;
        ray.intersection.tex_u = (atan2(n[0], n[2])/(2*M_PI)) + 0.5;
        ray.intersection.tex_v = n[1] * 0.5 + 0.5;
        return true;
    }

    return false;
}

bool intersects_top_or_bot(double t, Point3D origin, Vector3D direction,
                           Intersection intersection) {
    Point3D point = origin + t * direction;
    if (point[0] * point[0] + point[2] * point[2] < 1 && t > 0 && t < intersection.t_value) {
        return true;
    }
    return false;
}

bool UnitCylinder::intersect(Ray3D &ray, const Matrix4x4 &worldToModel, const Matrix4x4 &modelToWorld) {
    // Setup new ray origin and dir in model space
    Point3D ray_origin = worldToModel * ray.origin;
    Vector3D ray_dir = worldToModel * ray.dir;

    // Describe cylinder geometry here, center, zmin, zmax, radius
    Point3D center(0, 0, 0);
    double zmin = -0.5, zmax = 0.5;
    double radius = 1.0;

    // Setup values for solving quardratic
    double temp1, temp2, t, a, b, c, d, x1, x2;
    Vector3D normal;
    double err = 0.0001;

    Point3D intersectionPoint;

    if (ray_dir[2] != 0) {
        // Check if the ray intersects the two caps zmin and zmax
        temp1 = (zmin - ray_origin[2]) / ray_dir[2];
        temp2 = (zmax - ray_origin[2]) / ray_dir[2];

        // first set the intersect 
        t = temp2;
        Point3D _normal(0, 0, 1);

        if (temp1 < temp2) {
            t = temp1;
            // Construct the normal for the closer disk, which points on the
            // negative z axis
            Point3D _normal(0, 0, -1);
        }
        normal = _normal - center;
        normal.normalize();

        if (t <= err) {
            return false;
        }

        intersectionPoint = ray_origin + t * ray_dir;

        // Intersection with cap and bottom
        if (intersectionPoint[0] * intersectionPoint[0] + intersectionPoint[1] * intersectionPoint[1] <= radius * radius) {

            if (ray.intersection.none || t < ray.intersection.t_value) {
                ray.intersection.point = intersectionPoint;
                ray.intersection.normal = normal;
                ray.intersection.t_value = t;
                ray.intersection.none = false;

                return true;
            }
            return false;
        }

        // Intersection with body
        a = ray_dir[0] * ray_dir[0] + ray_dir[1] * ray_dir[1];
        b = ray_origin[0] * ray_dir[0] + ray_origin[1] * ray_dir[1];
        c = ray_origin[0] * ray_origin[0] + ray_origin[1] * ray_origin[1] - radius * radius;
        d = b * b - a * c;

        if (d < 0) {
            return false;
        }

        x1 = (-b + sqrt(d)) / a;
        x2 = (-b - sqrt(d)) / a;

        if (x1 > 0 || x2 > 0) {

            if (x1 > 0 && x2 > 0) {

                if (x1 <= x2) {
                    t = x1;
                }

                else {
                    t = x2;
                }
            }

            else if (x1 > 0) {
                t = x1;
            }

            else {
                t = x2;
            }
            intersectionPoint = ray_origin + t * ray_dir;
            normal[0] = intersectionPoint[0];
            normal[1] = intersectionPoint[1];
            normal[2] = 0;
            normal.normalize();

            if (intersectionPoint[2] < zmax && intersectionPoint[2] > zmin) {

                if (ray.intersection.none || t < ray.intersection.t_value) {
                    ray.intersection.point = modelToWorld * intersectionPoint;
                    Point3D normalTemp;
                    normalTemp[0] = intersectionPoint[0];
                    normalTemp[1] = intersectionPoint[1];
                    normalTemp[2] = 0;
                    ray.intersection.normal = modelToWorld * (normalTemp - center);
                    ray.intersection.t_value = t;
                    ray.intersection.none = false;
                    return true;
                }
                return false;
            }
            return false;
        }
        return false;
    }

    // ray_origin[2] == 0

    a = pow(ray_dir[0], 2) + pow(ray_dir[1], 2);
    b = 2 * (ray_dir[0] * ray_origin[0] + ray_dir[1] * ray_origin[1]);
    c = pow(ray_origin[0], 2) + pow(ray_origin[1], 2) - 1;
    d = b * b - 4 * a * c;

    if (d < 0) {
        return false;
    }

    x1 = (-b + sqrt(d)) / (2 * a);
    x2 = (-b - sqrt(d)) / (2 * a);

    if (x1 > 0 || x2 > 0) {

        if (x1 > 0 && x2 > 0) {

            if (x1 <= x2) {
                t = x1;
            }

            else {
                t = x2;
            }
        }

        else if (x1 > 0) {
            t = x1;
        }

        else {
            t = x2;
        }

        intersectionPoint = ray_origin + t * ray_dir;
        normal[0] = intersectionPoint[0];
        normal[1] = intersectionPoint[1];
        normal[2] = 0;
        normal.normalize();

        if (intersectionPoint[2] < zmax && intersectionPoint[2] > zmin) {

            if (ray.intersection.none || t < ray.intersection.t_value) {
                ray.intersection.point = modelToWorld * intersectionPoint;
                Point3D normalTemp;
                normalTemp[0] = intersectionPoint[0];
                normalTemp[1] = intersectionPoint[1];
                normalTemp[2] = 0;
                ray.intersection.normal = modelToWorld * (normalTemp - center);
                ray.intersection.t_value = t;
                ray.intersection.none = false;
                return true;
            }

            return false;
        }

        return false;
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
