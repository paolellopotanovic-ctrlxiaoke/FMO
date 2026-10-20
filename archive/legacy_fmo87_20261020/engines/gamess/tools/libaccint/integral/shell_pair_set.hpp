/*
 * shell_pair_set.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_INTEGRAL_SHELL_PAIR_SET_HPP_
#define LIBACCINT_INCLUDE_INTEGRAL_SHELL_PAIR_SET_HPP_
 
// External header(s)

// Internal header(s)
#include "core/shell_set.hpp"

// STL header(s)
#include <vector>
#include <tuple>

// Convenience typedef(s)
typedef libaccint::host::core::ShellSet ShellSet;

namespace libaccint {

    namespace host {

        namespace integral {

            /*! \class ShellPairSet
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief Set of Shell pairs
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            class ShellPairSet {

            public:

                /***** Constructor(s) and Destructors(s) *****/

                // Default constructor
                ShellPairSet();

                // Standard constructor
                ShellPairSet(ShellSet& shell_set_1, ShellSet& shell_set_2);

                // Copy constructor
                ShellPairSet(const ShellPairSet& shell_pair_set_to_copy);

                // Move constructor
                ShellPairSet(ShellPairSet&& shell_pair_set_to_move);

                // Default destructor
                ~ShellPairSet();

                /***** Public Member Functions *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Operator Overloads *****/

                // Copy assignment operator
                ShellPairSet& operator= (const ShellPairSet& shell_pair_set_to_copy);

                // Move assignment operator
                ShellPairSet& operator= (ShellPairSet&& shell_pair_set_to_move);

                // TODO: ADD_FOR_BETA_RELEASE

            private:

                /***** Private Member Functions *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Private Member Variables *****/

                // Angular momentum (i.e. shell type) of the 'a' shells
                int a_shell_angular_momenta_;

                // Angular momentum (i.e. shell type) of the 'b' shells
                int b_shell_angular_momenta_;

                // Degree of contraction for each of the 'a' shells
                int a_shells_num_primitives_;

                // Degree of contraction for each of the 'b' shells
                int b_shells_num_primitives_;

                // Vector of a,b shell index pairs
                std::vector< std::tuple<int, int> > ab_shell_pair_indices_;
                
                // TODO: Etc.

            };

        } // namespace libaccint::host::integral

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_INTEGRAL_SHELL_PAIR_SET_HPP_ */