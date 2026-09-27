/*
 * shell.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_CORE_SHELL_HPP_
#define LIBACCINT_INCLUDE_CORE_SHELL_HPP_

// External header(s)

// Internal header(s)

// STL header(s)
#include <array>

namespace libaccint {

    namespace host {

        namespace core {

            /*! \class Shell
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief Shell class
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            class Shell {

            public:

                /***** Constructor(s) and Destructors(s) *****/

                // Default constructor
                Shell();

                // Standard constructor
                Shell(/*TODO: PARAMS*/);

                // Copy constructor
                Shell(const Shell& shell_to_copy);

                // Move constructor
                Shell(Shell&& shell_to_move);

                // Default destructor
                ~Shell();

                /***** Public Member Functions *****/

                /***** Accessors *****/

                // Return the number of primitives in the Shell
                unsigned int num_primitives() const;

                // Return the angular momentum of the Shell
                unsigned int angular_momentum() const;

                // Return the expansion center of the Shell
                std::array<double, 3> center() const;

                // Return the set of contraction coefficients
                std::array<double, num_primitives> cont_coeffs() const;

                // Return the exponential factors
                std::array<double, num_primitives> exponents() const;

                // Return a particular contraction coefficient
                double cont_coeff(unsigned int cont_coeff_index) const;

                // Return a particular exponent
                double exponent(unsigned int exponent_index) const;

                /***** Setters *****/

                // Set a new expansion center for the Shell
                void center(std::array<double, 3> new_center);

                /***** Operator Overloads *****/

                // Copy assignment operator
                Shell& operator= (const Shell& shell_to_copy);

                // Move assignment operator
                Shell& operator= (Shell&& shell_to_move);

                // TODO: ADD_FOR_BETA_RELEASE

            private:

                /***** Private Member Functions *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Private Member Variables *****/

                // Angular momentum (i.e. shell type) of the Shell instance
                // (Note: Negative momentum values are reserved for hybrid shells (e.g. sp a.k.a. L shell))
                int angular_momentum_;

                // Expansion center of the shell
                std::array<double, 3> expansion_center_;

                // Array of exponential factors
                std::vector<double> exponents_;

                // Array of contraction coefficients
                std::vector<double> contraction_coefficients_;

            };

        } // namespace libaccint::host::core

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_CORE_SHELL_HPP_ */