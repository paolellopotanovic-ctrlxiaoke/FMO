/*
 * basis_set.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_CORE_BASIS_SET_HPP_
#define LIBACCINT_INCLUDE_CORE_BASIS_SET_HPP_

// External header(s)

// Internal header(s)
#include "atom.hpp"
#include "shell.hpp"
#include "shell_set.hpp"
#include "source/host/core/raw_basis_set.hpp"

// STL header(s)
#include <string>
#include <vector>

namespace libaccint {

    namespace host {

        namespace core {

            /*! \class BasisSet
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief The basis set class
            *
            *  \details TODO
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            class BasisSet {

            public:

                /***** Constructor(s) and Destructor(s) *****/

                  // Default constructor
                BasisSet();

                // Sorted BasisSet constructor
                BasisSet(std::string basis_name, std::vector<Atom>& atoms, bool sorted = true);

                // Copy constructor
                BasisSet(const BasisSet& basis_set_to_copy);

                // Move constructor
                BasisSet(BasisSet&& basis_set_to_move);

                // Default destructor
                ~BasisSet();

                /***** Public Member Functions *****/

                // Return the highest angular momentum value for a shell (and therefore a shell set)
                int max_ang_mom() const;

                // Return the number of shell sets
                unsigned int num_shell_sets() const;

                // Return the number of shells
                unsigned int num_shells() const;

                // Return the number of basis functions
                unsigned int num_basis_fns() const;

                /***** Operator Overloads *****/

                // Return the ShellSet with index 'shell_set_index' in the basis set
                // Note: A particular Shell instance may be accessed by basis_set[shell_set_index][shell_index], since
                //       ShellSet has its own operator[] overload which returns a Shell instance
                ShellSet& operator[] (unsigned int shell_set_index);

                // Const version of ShellSet accessor
                const ShellSet& operator[] (unsigned int shell_set_index) const;

                // Copy assignment operator
                BasisSet& operator= (const BasisSet& basis_set_to_copy);

                // Move assignment operator
                BasisSet& operator= (BasisSet&& basis_set_to_move);

                // TODO: Add further overloads as necessary

            private:

                /***** Private Member Functions *****/

                  // TODO: ADD_FOR_BETA

                /***** Private Member Variables *****/

                // Underlying raw basis set data
                RawBasisSet& raw_basis_set_;

            };

        } // namespace libaccint::host::core

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_CORE_BASIS_SET_HPP_ */