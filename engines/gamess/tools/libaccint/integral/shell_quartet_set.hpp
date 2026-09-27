/*
 * shell_quartet_set.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_INTEGRAL_SHELL_QUARTET_SET_HPP_
#define LIBACCINT_INCLUDE_INTEGRAL_SHELL_QUARTET_SET_HPP_

// External header(s)

// Internal header(s)
#include "shell_pair_set.hpp"

// STL header(s)
#include <vector>
#include <array>

namespace libaccint {

    namespace host {

        namespace integral {

            /*! \class ShellQuartetSet
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief Set of Shell quartets. Used when screening is employed
            *         on the host-side.
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            class ShellQuartetSet {

            public:

                /***** Constructor(s) and Destructors(s) *****/

                // Default constructor
                ShellQuartetSet();

                // Standard constructor. Performs screening to reduce the full number of shell quartets
                // in the provided arguments to only those which are significant.
                ShellQuartetSet(ShellPairSet& shell_pair_set_1, ShellPairSet& shell_pair_set_2);

                // Copy constructor
                ShellQuartetSet(const ShellQuartetSet& shell_quartet_set_to_copy);

                // Move constructor
                ShellQuartetSet(ShellQuartetSet&& shell_quartet_set_to_move);

                // Default destructor
                ~ShellQuartetSet();

                /***** Public Member Functions *****/

                // Get the number of shell quartets in the set
                unsigned int num_quartets() const;

                // Get the set of angular momentum values for the shells in the set
                // (e.g. a set of shells for the class (ds|pp) would return {2, 0, 1, 1})
                std::array<int, 4> ang_mom_values() const;

                // Gett the degree of contraction for each of the shell types in the quartets
                std::array<unsigned int, 4> degrees_of_cont() const;

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Operator Overloads *****/

                // Obtain a given quartet of shell indices
                std::array<int, 4> operator[] (unsigned int shell_quartet_index);

                // Copy assignment operator
                ShellQuartetSet& operator= (const ShellQuartetSet& shell_quartet_set_to_copy);

                // Move assignment operator
                ShellQuartetSet& operator= (ShellQuartetSet&& shell_quartet_set_to_move);

                // TODO: ADD_FOR_BETA_RELEASE

            private:

                /***** Private Member Functions *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Private Member Variables *****/

                // Angular momentum (i.e. shell type) of the 'a' shells
                int a_shell_angular_momenta_;

                // Angular momentum (i.e. shell type) of the 'b' shells
                int b_shell_angular_momenta_;

                // Angular momentum (i.e. shell type) of the 'c' shells
                int c_shell_angular_momenta_;

                // Angular momentum (i.e. shell type) of the 'd' shells
                int d_shell_angular_momenta_;

                // Degree of contraction for each of the 'a' shells
                int a_shells_num_primitives_;

                // Degree of contraction for each of the 'b' shells
                int b_shells_num_primitives_;

                // Degree of contraction for each of the 'c' shells
                int c_shells_num_primitives_;

                // Degree of contraction for each of the 'd' shells
                int d_shells_num_primitives_;

                // Vector of a,b shell index pairs
                std::vector< std::array<int, 4> > abcd_shell_quartet_indices_;

                // TODO: ADD_FOR_BETA_RELEASE

            };

        } // namespace libaccint::host::integral

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_INTEGRAL_SHELL_QUARTET_SET_HPP_ */