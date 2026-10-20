/*
 * shell_set.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_CORE_SHELL_SET_HPP_
#define LIBACCINT_INCLUDE_CORE_SHELL_SET_HPP_

// External header(s)

// Internal header(s)
#include "shell.hpp"

// STL header(s)
#include <array>

namespace libaccint {

    namespace host {

        namespace core {

            /*! \class ShellSet
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief ShellSet class
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            class ShellSet {

            public:

                /***** Constructor(s) and Destructors(s) *****/

                // Default constructor
                ShellSet();

                // Standard constructor
                ShellSet(/*TODO: PARAMS*/);

                // Copy constructor
                ShellSet(const ShellSet& shell_set_to_copy);

                // Move constructor
                ShellSet(ShellSet&& shell_set_to_move);

                // Default destructor
                ~ShellSet();

                /***** Public Member Functions *****/

                // Number of shells within the shell set
                unsigned int num_shells() const;

                // Return the angular momentum of all shells in the shell set
                int angular_momentum() const;

                // Return the number of primitives for each shell in the shell set
                unsigned int num_primitives() const;

                /***** Operator Overloads *****/

                // Shell accessor
                Shell& operator[] (unsigned int shell_index);

                // Const version of Shell accessor
                const Shell& operator[] (unsigned int shell_index) const;

                // Copy assignment operator
                ShellSet& operator= (const ShellSet& shell_set_to_copy);

                // Move assignment operator
                ShellSet& operator= (ShellSet&& shell_set_to_move);

                // TODO: ADD_FOR_BETA_RELEASE

            private:

                /***** Private Member Functions *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Private Member Variables *****/

                // Angular momentum (i.e. shell type) of all Shell instances in the set
                const unsigned int angular_momentum_;

                // Number of primitives of all Shell instances in the set
                const unsigned int num_primitives_;

                // Vector of Shell instances
                std::vector<Shell> shells_;

                // TODO: ADD_FOR_BETA_RELEASE

            };

        } // namespace libaccint::host::core

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_CORE_SHELL_SET_HPP_ */