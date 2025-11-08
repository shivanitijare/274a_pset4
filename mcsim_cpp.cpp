#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <chrono> // for generating the seed
#include <random> // for random numbers

typedef std::array<double, 3> AtomCoord;
typedef std::vector<AtomCoord> Coordinates;

// A Global! Probably shouldn't be used in real code
std::default_random_engine re;

/*! Generate a random double within a given range */
double random_double(double lower_bound, double upper_bound)
{
    std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
    return dist(re);
}


/*! Generate a random integer within a given range
    The generated integer will be on the range [a,b)
*/
double random_integer(int lower_bound, int upper_bound)
{
    // dist will return [a,b] but we want [a,b)
    std::uniform_int_distribution<int> dist(lower_bound, upper_bound - 1);
    return dist(re);
}


double calculate_LJ(double r_ij)
{
    double inv_r_ij = 1.0 / r_ij;
    double r6_term = pow(inv_r_ij, 6.0);
    double r12_term = r6_term * r6_term;
    double pairwise_energy = 4.0 * (r12_term - r6_term);
    return pairwise_energy;
}


double calculate_distance(AtomCoord coord1, AtomCoord coord2, double box_length = -1.0)
{

    double distance = 0;

    for(int i = 0; i < 3; i++) {
        double dim_dist = (coord1[i] - coord2[i]);

        if (box_length != -1.0) {
            dim_dist = dim_dist - box_length * round(dim_dist / box_length);
        }
        dim_dist = dim_dist * dim_dist;
        distance += dim_dist;
    }

    distance = sqrt(distance);

    return distance;
}


double calculate_total_energy(Coordinates coordinates, double box_length, double cutoff)
{

    double total_energy = 0.0;

    for(int i = 0; i < coordinates.size(); i++)
    {
        for(int j = i + 1; j < coordinates.size(); j++)
        {

            double dist_ij = calculate_distance(coordinates[i], coordinates[j], box_length);

            if (dist_ij < cutoff)
                total_energy += calculate_LJ(dist_ij);
        }
    }

    return total_energy;
}


double calculate_tail_correction(int num_particles, double box_length, double cutoff)
{
    double const1 = (8.0 * M_PI * num_particles * num_particles) / (3 * pow(box_length, 3));
    double const2 = (1.0 / 3.0) * pow((1.0 / cutoff), 9) - pow((1.0 / cutoff), 3);
    return const1 * const2;
}


double calculate_pair_energy(Coordinates coordinates, int i_particle, double box_length, double cutoff)
{

    double e_total = 0.0;
    int num_atoms = coordinates.size();

    for(int j_particle = 0; j_particle < num_atoms; j_particle++)
    {
        if (j_particle != i_particle)
        {
            double dist_ij = calculate_distance(coordinates[i_particle], coordinates[j_particle], box_length);

            if (dist_ij < cutoff)
            {
                double interaction_energy = calculate_LJ(dist_ij);
                e_total += interaction_energy;
            }
        }
    }

    return e_total;
}


bool accept_or_reject(double delta_e, double beta)
{

    bool accept = false;
    if (delta_e <= 0.0)
        accept = true;
    else
    {
        double random_number = random_double(0, 1);
        double p_acc = exp(-beta * delta_e);
        if (random_number < p_acc)
            accept = true;
        else
            accept = false;
    }

    return accept;
}


std::pair<Coordinates, double> read_xyz(std::string file_path)
{
    // Opens up a file stream forinput
    std::ifstream infile(file_path);

    // Check that it was successfully opened
    if (!infile.is_open()) {
        throw std::runtime_error("File path in read_xyz does not exist!");
    }

    double dummy; // Data that is thrown away (box length, atom indices)
    double box_length;
    int num_atoms;

    // Grab box_length from first number, throw the rest away
    infile >> box_length >> dummy >> dummy;

    // now the number of atoms
    infile >> num_atoms;

    // Uncomment to help troubleshoot
    // std::cout << "Box length: " << box_length << " natoms: " << num_atoms << std::endl;

    // Stores the atomic coordinates
    // Remember, this is a vector of arrays
    Coordinates coords;

    for(int i = 0; i < num_atoms; i++)
    {
        AtomCoord coord;

        // Throws away the atom index
        infile >> dummy >> coord[0] >> coord[1] >> coord[2];

        // Add to the vector
        coords.push_back(coord);
    }

    // Makes an appropriate pair object
    return std::make_pair(coords, box_length);
}


Coordinates run_simulation(
    Coordinates coordinates,
    const double box_length,
    const double cutoff,
    const double reduced_temperature,
    const int num_steps,
    const double max_displacement = 0.1,
    const int freq = 1000)
{

    std::vector<int> steps;
    std::vector<double> energies;
    std::vector<Coordinates> all_coordinates;

    double beta = 1.0 / reduced_temperature;

    int num_particles = coordinates.size();

    double total_energy = calculate_total_energy(coordinates, box_length, cutoff);
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff);

    for(int i = 0; i < num_steps; i++)
    {

        int random_particle = random_integer(0, num_particles);

        double current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double x_rand = random_double(-max_displacement, max_displacement);
        double y_rand = random_double(-max_displacement, max_displacement);
        double z_rand = random_double(-max_displacement, max_displacement);

        coordinates[random_particle][0] += x_rand;
        coordinates[random_particle][1] += y_rand;
        coordinates[random_particle][2] += z_rand;

        double proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double delta_energy = proposed_energy - current_energy;

        bool accept = accept_or_reject(delta_energy, beta);

        if (accept)
            total_energy += delta_energy;
        else
        {
            coordinates[random_particle][0] -= x_rand;
            coordinates[random_particle][1] -= y_rand;
            coordinates[random_particle][2] -= z_rand;
        }

        if (i % freq == 0)
        {
            std::cout << i << ' ' << total_energy / num_particles << std::endl;
            steps.push_back(i);
            energies.push_back(total_energy / num_particles);
            all_coordinates.push_back(coordinates);
        }
    }

    return coordinates;
}


int main(void)
{
    double reduced_temperature = 1.5;
    double max_displacement = 0.1; 
    double cutoff = 3.0;

    int num_steps = 50000;

    re.seed(std::chrono::system_clock::now().time_since_epoch().count());

    std::pair<Coordinates, double> xyz_info = read_xyz("lj_sample_config_periodic1.txt");
    Coordinates coords = xyz_info.first;
    double box_length = xyz_info.second;

    Coordinates newcoords = run_simulation(coords, box_length, cutoff, reduced_temperature, num_steps);

    return 0;
}
