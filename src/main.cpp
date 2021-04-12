#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <random>
#include <limits>
#include <memory> // unique_ptr

#include "tgaimage.h"
#include "geometry.h"

constexpr int      g_depth_cutoff = 8;
constexpr int      g_screen_x     = 1024;
constexpr int      g_screen_y     = 768;
constexpr double   g_view_ang     = 0.5135;
const Vector3d g_cam_pos      = Vector3d(50, 52, 295.6);
const Vector3d g_cam_dir      = Vector3d(0, -0.042612, -1).normalize();
// the number 140 pushes camera rays forward so they skip any intersections
// close to the camera
constexpr int    g_ray_push     = 140;

std::random_device seed;
std::default_random_engine generator(seed());
std::uniform_real_distribution<double> distr(0.0,1.0);

struct Ray {
    Vector3d origin;
    Vector3d direction;
    Ray(Vector3d o_, Vector3d d_) : origin(o_), direction(d_) {};
};

// Material types
enum Refl_t {
    DIFFUSE,
    SPECULAR,
    REFRACTIVE,
};

class Shape {
	public: 
		Vector3d position;
		Vector3d emission;
		Vector3d color;
		Refl_t refl;
		virtual double intersect(const Ray &ray) const = 0;
		virtual Vector3d normal_v(const Vector3d& intersect) const = 0;
		virtual Vector3d get_pos() const = 0;
		virtual Vector3d get_emission() const = 0;
		virtual Vector3d get_color(const Vector3d& intersect_v) const = 0;
		virtual Refl_t get_refl() const = 0;
        virtual ~Shape() = default;
};

class Sphere : public Shape {
	private:
        double radius;
    public:
        Vector3d position;
        Vector3d emission;
        Vector3d color;
        Refl_t refl; // Reflection type
        
        double intersect(const Ray& ray) const override;
		Vector3d normal_v(const Vector3d& intersect) const override;
		Vector3d get_pos() const override {return position;};
		Vector3d get_emission() const override {return emission;};
		Vector3d get_color(const Vector3d& intersect_v) const override;
		Refl_t get_refl() const override {return refl;};
		void load_texture(const std::string& filename, TGAImage &image);
		std::string texturefile;
		TGAImage diffusemap_;
		TGAColor diffuse_lookup(double u, double v) const {
			return diffusemap_.get(static_cast<int>(u * diffusemap_.get_width()), static_cast<int>(v * diffusemap_.get_height()));
		}
		Sphere(double radius_, Vector3d position_, Vector3d emission_, Vector3d color_, Refl_t refl_, std::string texturefile_="NONE") :
               		radius(radius_), position(position_), emission(emission_), color(color_), refl(refl_), texturefile(texturefile_){
						if (texturefile != "NONE") {
							load_texture(texturefile, diffusemap_);
						}
					};
        ~Sphere() override = default;
};

void Sphere::load_texture(const std::string& filename, TGAImage& image){
	std::cout << "Loading texture " << filename << "\n";
	image.read_tga_file(filename.c_str());
	// image.flip_vertically();
	texturefile = filename;
}

std::ostream& operator<<(std::ostream& out, const Vector3d& v) {
    out << "[" << v.x << ", " << v.y << ", " << v.z << "]";
    return out;
}
std::ostream& operator<<(std::ostream& out, const Sphere& sph) {
    out << "p:" << sph.position << ", col: " << sph.color << ", em: " << sph.emission << ", refl: " << sph.refl << "\n";
    return out;
}

std::ostream& operator<<(std::ostream& out, const Ray& r) {
    out << r.origin << " + L" << r.direction;
    return out;
}

Vector3d Sphere::get_color(const Vector3d& intersect_v) const {
	if (texturefile != "NONE") {
		Vector3d vec = (intersect_v - position).normalize();
		double u = 0.5 + std::atan2(vec.x, vec.z) / (2 * M_PI);
		double v = 0.5 - std::asin(vec.y) / M_PI;
	
		TGAColor texture_color = diffuse_lookup(u, v);
		return Vector3d(texture_color.bgra[2] / 255.0, texture_color.bgra[1] / 255.0, texture_color.bgra[0] / 255.0);
	}
    
    return color;
}

inline double clamp(double x) {
	// Clamps a value to between 0 and 1
    if (x < 0) {
        return 0;
    }
    if (x > 1) {
        return 1;
    }
    return x;
}

int gamma_corr(double x) {
    return static_cast<int>(pow(clamp(x), 1/2.2) * 255 + 0.5);
}

double Sphere::intersect(const Ray &ray) const{
    // Tests if this sphere intersects with the specified ray
    // Returns 0 if there is no hit,
    // and returns the distance to the closest pont of intersection
    
    // Solving the equation 
    // distance ** 2 (dir.dir) + 2 * distance.(origin_to_point) + (origin_to_point).(origin_to_point) = radius**2
    Vector3d origin_to_pos = position - ray.origin;
    double epsilon = 1e-4;
    double b = origin_to_pos.dot(ray.direction);
    double det = (b * b) - origin_to_pos.dot(origin_to_pos) + (radius * radius);
    if (det < 0) {
        return -1;
    }
    
    det = std::sqrt(det);
  
    if (b - det > 0 && (b - det) > epsilon) {
		return b - det;
	}
    
    if (b + det > 0 && (b + det) > epsilon) {
		return b + det;
	}
    
    return -1;
}

Vector3d Sphere::normal_v(const Vector3d& intersect) const{
	return (intersect - position).normalize();
}

std::vector<std::unique_ptr<Shape>> load_scene_file(const std::string& filename) {
	// Reads in a vector of spheres from a specified
	// input file.
	std::ifstream infile(filename);
	std::vector<std::unique_ptr<Shape>> scene;
	if (infile.is_open()) {
		std::string line;
		while ( std::getline(infile, line) ) {
			if (line.rfind("//", 0) == 0) {
				continue;
			}
			double radius;
			Vector3d pos;
			Vector3d emission;
			Vector3d color;
			std::string refl_t_str;
			Refl_t refl_type;
			char delim;
			
			std::istringstream iss(line);
			iss >> radius >> delim;
			iss >> pos.x >> pos.y >> pos.z >> delim;
			iss >> emission.x >> emission.y >> emission.z >> delim;
			iss >> color.x >> color.y >> color.z >> delim;
			iss >> refl_t_str >> delim;
			if (refl_t_str == "DIFFUSE") {
				refl_type = DIFFUSE;
			} else if (refl_t_str == "REFRACTIVE") {
				refl_type = REFRACTIVE;				
			} else if (refl_t_str == "SPECULAR") {
				refl_type = SPECULAR;				
			} else {
				refl_type = DIFFUSE;
				std::cerr << "Bad reflection type: " << refl_t_str << ". Using default of DIFFUSE.\n";
				
			}
			std::string texture_file;
			iss >> texture_file;
			std::cout << texture_file << "\n";
			scene.emplace_back(std::make_unique<Sphere>(Sphere(radius, pos, emission, color, refl_type, texture_file)));
		}
		infile.close();		
	} else {
		std::cerr << "Could not read " << filename << "\n";
	}
	return scene;
}

Vector3d get_diffuse_color(const Vector3d& normal_v, const Vector3d& normal_v_oriented, const Vector3d& intersect_v, int depth, const std::vector<std::unique_ptr<Shape>>& scene);
Vector3d get_specular_color(const Ray& ray, const Vector3d& normal_v, const Vector3d& intersect_v, int depth, const std::vector<std::unique_ptr<Shape>>& scene);

Vector3d radiance(const Ray& ray, int depth, const std::vector<std::unique_ptr<Shape>>& scene) {
	// Recursively calculates the colour of the ray we've sent out.
	// Rays can have three types of intersection:
	// - Diffuse, scattered in a random hemisphere oriented by the normal vector
	// - Specular, with angle of incidence = angle of reflection
	// - Refractive, which splits into a reflected and transmitted ray
    double shortest_dist = 1e20; // Distance to intersection
    int id = 0; // Index of the intersected object
    
    for (uint32_t i = 0; i < scene.size(); ++i) {
        double dist = scene[i]->intersect(ray);
        if (dist < shortest_dist && dist > 0) {
            shortest_dist = dist;
            id = i;
        }
    }
    
    if (shortest_dist == 1e20) {
        // If there's nothing to hit, return black
        return Vector3d(0, 0, 0);
    }
    
    // Increment the reflection depth, 
	// and test our depth
    depth++;
    if (depth > g_depth_cutoff) {
        return scene[id]->get_emission();
    }  
    // Set up the vectors describing the intersection
    Vector3d intersect_v = ray.origin + (ray.direction * shortest_dist);
    Vector3d normal_v = scene[id]->normal_v(intersect_v);
    Vector3d normal_v_oriented;
    if (normal_v.dot(ray.direction) < 0) {
        normal_v_oriented = normal_v;
    } else {
		normal_v_oriented = normal_v * -1;
	}
    
    Vector3d color = scene[id]->get_color(intersect_v);

    switch (scene[id]->get_refl()) {
        case DIFFUSE:
        {							
			Vector3d next_color = get_diffuse_color(normal_v, normal_v_oriented, intersect_v, depth, scene);
            return scene[id]->get_emission() + mult(color, next_color);
        }
        case SPECULAR:
        {
			Vector3d next_color = get_specular_color(ray, normal_v, intersect_v, depth, scene);
            return scene[id]->get_emission() + mult(color, next_color);
        }
        case REFRACTIVE:
        {
            Vector3d new_dir = ray.direction - normal_v * 2 * normal_v.dot(ray.direction);
            Ray refl_ray(intersect_v, new_dir);
			
			// Is this a ray from outside the sphere going in, or not?
            bool into = normal_v.dot(normal_v_oriented) > 0;
        
			// Indices of refaction
            double ior_air = 1.0;
            double ior_glass = 1.5;
			double nnt;
			if (into) {
				nnt = ior_air / ior_glass;
			} else {
				nnt = ior_glass / ior_air;
			}
            double ddn = normal_v_oriented.dot(ray.direction);
        
            double cos2t = 1 - ((nnt * nnt) * (1 - ddn * ddn));
            if (cos2t < 0) {
                // Total internal reflection case
                return scene[id]->get_emission() + mult(color, radiance(Ray(intersect_v, new_dir), depth, scene));
            }

			// R0 is the reflection at "Normal Incidence"
            double R0 = ((nnt - 1) * (nnt - 1)) / ((nnt + 1) * (nnt + 1));
            double c;
            Vector3d tdir;
            if (into) {
                tdir = (ray.direction * nnt + normal_v * (ddn * nnt + sqrt(cos2t))).normalize(); 
                c = 1 + ddn;
            } else {
                tdir = (ray.direction * nnt - normal_v * (ddn * nnt + sqrt(cos2t))).normalize();
                c = 1 - tdir.dot(normal_v);
            }
        
			// But at other angles, it is given by
			// R0 + (1 - R0)(1 - cos(theta))^5
            double Re = R0 + ((1 - R0) * std::pow(c, 5.0));
			// Transmitted irradiance via conservation of energy
            double Tr = 1 - Re;
            double P  = 0.25 + (0.5 * Re);
            double RP = Re / P;
            double TP = Tr / (1 - P);
        
            if (distr(generator) < P) {
                return scene[id]->get_emission() + mult(color, radiance(refl_ray, depth, scene)) * RP;
            }
            return scene[id]->get_emission() + mult(color, radiance(Ray(intersect_v, tdir), depth, scene)) * TP;
        }
        default:
            std::cerr << "ERROR - END OF SWITCH BLOCK\n";
            return Vector3d(-1, -1, -1);
    }

}

Vector3d get_diffuse_color(const Vector3d& normal_v, const Vector3d& normal_v_oriented, const Vector3d& intersect_v, int depth, const std::vector<std::unique_ptr<Shape>>& scene) {
    // Ideal diffuse reflection, with
    // totally random bounces, sampling a hemisphere
    // at the point of incidence
	
    double r1 = 2 * M_PI * distr(generator);
    double r2 = distr(generator);
    Vector3d u;
        
    if (std::abs(normal_v.x) > 0.1) {
        u = cross(Vector3d(0, 1, 0), normal_v_oriented).normalize();
    } else {
         u = cross(Vector3d(1, 0, 0), normal_v_oriented).normalize();           
    }
        
    Vector3d v = cross(normal_v_oriented, u);
    Vector3d new_dir = (u * std::cos(r1) * std::sqrt(r2) +
                        v * std::sin(r1) * std::sqrt(r2) +
                        normal_v_oriented * std::sqrt(1 - r2)).normalize();
								
	return radiance(Ray(intersect_v, new_dir), depth, scene);
}

Vector3d get_specular_color(const Ray& ray, const Vector3d& normal_v, const Vector3d& intersect_v, int depth, const std::vector<std::unique_ptr<Shape>>& scene){
            // Ideal specular reflection, with 
            // angle of incidence = angle of reflection.
            // Add the emissive color of this sphere along with
            // whatever colours we pick up from the bounces
	Vector3d new_dir = ray.direction - normal_v * (2 * normal_v.dot(ray.direction));
    return radiance(Ray(intersect_v, new_dir), depth, scene);
}
int main(int argc, char** argv) {

    int samples = 64;
	std::vector<std::unique_ptr<Shape>> scene;
    if (argc == 2) {
        std::string filename(argv[1]);
		std::cout << "Loading from " << filename << "\n";
		scene = load_scene_file(argv[1]);
	} else if (argc == 3) {
        std::string filename(argv[1]);
		std::cout << "Loading from " << filename << "\n";
		scene = load_scene_file(argv[1]);
		samples = atoi(argv[2]);
	} else {
		std::cerr << "Warning, no scene descriptor specified. Exiting.\n";
		std::exit(EXIT_FAILURE);
	}
    
    Ray camera(g_cam_pos, g_cam_dir);
    TGAImage image(g_screen_x, g_screen_y, TGAImage::RGB);
    Vector3d cam_x(g_screen_x * g_view_ang / g_screen_y, 0, 0);
    Vector3d cam_y = cross(cam_x, camera.direction).normalize() * g_view_ang;
    Vector3d rad(0, 0, 0);
	auto start_time = std::chrono::high_resolution_clock::now();
	int last_percent_done = 0;
//#pragma omp parallel for schedule(dynamic, 1) private(rad)       // OpenMP
    for (int y = 0; y < g_screen_y; ++y) {
        for (int x = 0; x < g_screen_x; ++x) {
            //std::cout << "Rendering pixel (" << x << ", " << y << ")\n";
            rad = Vector3d(0, 0, 0);  
            for (int s = 0; s < samples; ++s) {
                double r1 = 2 * distr(generator);
                double r2 = 2 * distr(generator);
                double dx;
                double dy;
                if (r1 < 1) {
                    dx = std::sqrt(r1) - 1;
                } else {
                    dx = 1 - sqrt(2 - r1);
                }
                    
                if (r2 < 1) {
                    dy = std::sqrt(r1) - 1;
                } else {
                    dy = 1 - std::sqrt(2 - r2);
                }
                    
                Vector3d dir = cam_x * (((dx / 2) + x)/g_screen_x - 0.5) + 
                               cam_y * (((dy / 2) + y)/g_screen_y - 0.5) +
                               camera.direction;
                dir = dir.normalize();

                rad = rad + radiance(Ray(camera.origin + (dir*g_ray_push), dir), 0, scene) / samples;
            }
            TGAColor color(gamma_corr(rad.x), gamma_corr(rad.y), gamma_corr(rad.z));
            image.set(x, y, color);
        }
		auto now_time = std::chrono::high_resolution_clock::now();
		int percent_done = 100 * y / g_screen_y;
		if (percent_done > last_percent_done) {
			last_percent_done = percent_done;
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(now_time - start_time);
			auto time_left = time_span * (g_screen_y - y) / y;
			
			auto int_time = static_cast<int>(time_left.count());
			int hours = int_time / 3600;
			int minutes = (int_time % 3600) / 60;
			int seconds = int_time % 60;
			if (hours > 0) {
				std::cout << "\r " << percent_done << "% complete. ETA: " << hours << "h" << minutes << "m" << seconds << "s.    ";	
			} else if (minutes > 0) {
				std::cout << "\r " << percent_done << "% complete. ETA: " << minutes << "m" << seconds << "s                     ";				
			} else {
				std::cout << "\r " << percent_done << "% complete. ETA: " << seconds << "s.                                      ";				
			}
		}
		

    }
    image.flip_vertically();
    image.write_tga_file("raytraced.tga");
	std::cout << "\nWritten to file " << "raytraced.tga";
	auto now_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(now_time - start_time);
	std::cout << " in " << time_span.count() << " seconds\n";
    return 0;
}

