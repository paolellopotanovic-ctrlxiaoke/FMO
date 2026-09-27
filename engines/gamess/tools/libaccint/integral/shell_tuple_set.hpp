/*
 * shell_tuple_set.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_INTEGRAL_SHELL_TUPLE_SET_HPP_
#define LIBACCINT_INCLUDE_INTEGRAL_SHELL_TUPLE_SET_HPP_

// External header(s)

// Internal header(s)
#include "core/shell_set.hpp"

// STL header(s)
#include <vector>
#include <array>

// Convenience typedef(s)
typedef libaccint::host::core::ShellSet ShellSet;

namespace libaccint {

    namespace host {

        namespace integral {

            /*! \class ShellTupleSet
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief Set of Shell tuples (e.g. pairs, quartets, etc.)
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            template<unsigned int tuple_size>
            class ShellTupleSet {

            public:

                /***** Constructor(s) and Destructors(s) *****/

                // Default constructor
                ShellTupleSet();

                // Standard constructor
                ShellTupleSet(std::vector<ShellSet>& shell_sets);

                // Copy constructor
                ShellTupleSet(const ShellTupleSet& shell_tuple_set_to_copy);

                // Move constructor
                ShellTupleSet(ShellTupleSet&& shell_tuple_set_to_move);

                // Default destructor
                ~ShellTupleSet();

                /***** Public Member Functions *****/

                // Get the number of shell tuples in the set
                unsigned int tuple_size() const;

                // Get the set of angular momentum values for the shells in the set
                std::array<int, tuple_size> ang_mom_values() const;

                // Get the degree of contraction for each of the shell types
                std::array<unsigned int, tuple_size> degrees_of_cont() const;

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Operator Overloads *****/

                // Obtain a given tuple of shell indices
                std::array<int, tuple_size> operator[] (unsigned int shell_tuple_index);

                // Copy assignment operator
                ShellTupleSet& operator= (const ShellTupleSet& shell_tuple_set_to_copy);

                // Move assignment operator
                ShellTupleSet& operator= (ShellTupleSet&& shell_tuple_set_to_move);

                // TODO: ADD_FOR_BETA_RELEASE

            private:

                /***** Private Member Functions *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Private Member Variables *****/

                // Angular momentum (i.e. shell type) of each shell set
                std::array<int, tuple_size> shell_set_angular_momenta_;

                // Degree of contraction of each shell set
                std::array<unsigned int, tuple_size> shell_set_degrees_of_contraction_;

                // Vector of shell set tuple indices
                std::vector< std::array<unsigned int, tuple_size> > shell_set_tuple_indices_;

                // TODO: Etc.

            };

        } // namespace libaccint::host::integral

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_INTEGRAL_SHELL_TUPLE_SET_HPP_ */