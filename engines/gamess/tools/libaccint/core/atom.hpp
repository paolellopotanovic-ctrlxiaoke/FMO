/*
 * atom.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_CORE_ATOM_HPP_
#define LIBACCINT_INCLUDE_CORE_ATOM_HPP_

// External header(s)

// Internal header(s)

// STL header(s)
#include <array>

namespace libaccint {

    namespace host {

        namespace core {

            /*! \struct Atom
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief basic Atom struct
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            struct Atom {

                /***** Constructor(s) and Destructor *****/

                // Default constructor
                Atom();

                // Standard constructor
                Atom(unsigned int atomic_number, std::array<double, 3> center) :
                    atomic_number_{ atomic_number },
                    center_{ center }
                {}

                // Copy constructor
                Atom(const Atom& atom_to_copy);

                // Move constructor
                Atom(Atom&& atom_to_move);

                // Default destructor
                ~Atom();

                /***** Public Member Functions *****/

                /***** Getters *****/

                // Return atomic number
                inline unsigned int atomic_number() const { return atomic_number_; }

                // Return center
                inline std::array<double, 3> center() const { return center_; }

                /***** Setters *****/

                // Update center
                void center(std::array<double, 3> new_center) { center_ = new_center; }

                /***** Operator Overloads *****/

                // Copy assignment operator
                Atom& operator= (const Atom& atom_to_copy);

                // Move assignment operator
                Atom& operator= (Atom&& atom_to_move);

                /***** Public Member Variables *****/

                // Atomic number
                const int atomic_number_;

                // Cartesian position
                std::array<double, 3> center_;

                // TODO: Add any other parameters as necessary

            };

        } // namespace libaccint::host::core

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_CORE_ATOM_HPP_ */
