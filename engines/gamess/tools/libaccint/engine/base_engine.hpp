/*
 * base_engine.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_ENGINE_BASE_ENGINE_HPP_
#define LIBACCINT_INCLUDE_ENGINE_BASE_ENGINE_HPP_

// External header(s)

// Internal header(s)
#include "core/basis_set.hpp"
#include "integral_operator.hpp"
#include "integral_consumer.hpp"

// STL header(s)
#include <vector>

// Convenience typedef(s)
typedef libaccint::host::core::BasisSet BasisSet;
typedef libaccint::host::integral::IntegralOperator IntegralOperator;

namespace libaccint {

    namespace host {

        namespace engine {

            /*! \class BaseEngine
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief Abstract base Engine class for common functionality/parameters
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            class BaseEngine {

            public:

                /***** Constructor(s) and Destructor(s) *****/

                // Default constructor
                BaseEngine();

                // Standard constructor
                BaseEngine(BasisSet& basis_set);

                // Standard constructor with specified integral operator
                BaseEngine(BasisSet& basis_set, IntegralOperator integral_operator);

                // Copy constructor
                BaseEngine(const BaseEngine& base_engine);

                // Move constructor
                BaseEngine(BaseEngine&& base_engine);

                // TODO: Etc.

                // Default destructor
                ~BaseEngine();

                // TODO: Etc.

                /***** Public Member Functions *****/

                /***** Accessors *****/
                // Get the basis set name
                std::string name() const;

                // Get the current integral operator
                IntegralOperator integral_operator() const;

                // Get the number of basis sets currently attached to the engine
                unsigned int num_basis_sets() const;

                // Detach a basis set from the engine
                void remove_basis_set(BasisSet* basis_set_to_remove);

                // Attach a new basis set to the engine
                void add_basis_set(BasisSet* basis_set_to_add);

                // Assign an integral operator
                void integral_operator(IntegralType new_operator);
                
                // TODO: Etc.

                /***** Integral Compute Driver Functions *****/

                // The compute function for integrals over shell pairs
                double* compute(unsigned int a_shell_index, unsigned int b_shell_index);

                // The compute function for integrals over shell triplets (Reserved for RI-type computations)
                double* compute(unsigned int a_shell_index, unsigned int b_shell_index, unsigned int c_shell_index);

                // The integral compute function for integrals over shell quartets
                double* compute(unsigned int a_shell_index, unsigned int b_shell_index,
                                unsigned int c_shell_index, unsigned int d_shell_index);

                // The compute function for integrals over ShellSet pairs
                double* compute(unsigned int a_shell_set_index, unsigned int b_shell_set_index);

                // The compute function for integrals over ShellSet triplets (Reserved for RI-type computations)
                double* compute(unsigned int a_shell_set_index, unsigned int b_shell_set_index, unsigned int c_shell_set_index);

                // The compute function for integrals over ShellSet quartets
                double* compute(unsigned int a_shell_set_index, unsigned int b_shell_set_index,
                                unsigned int c_shell_set_index, unsigned int d_shell_set_index);

                // Compute a set of integrals over a set of unique shell pairs (Used for screened shell pair sets)
                double* compute(ShellPairSet& shell_pair_set);

                // Compute a set of integrals over a set of unique shell quartets (Used for screened shell quartet sets)
                double* compute(ShellQuartetSet& shell_quartet_set);

                // Compute a set of integrals over a set of unique shell tuples (Used for screened shell tuple sets)
                double* compute(ShellTupleSet& shell_tuple_set);

                // TODO: Add functionality for computing integrals over Shells/ShellSets/Shell(Pair/Quartet/Tuple)Sets in different basis sets

                /***** Operator Overloads *****/

                // Copy assignment operator overload
                BaseEngine& operator=(const BaseEngine& base_engine);

                // Move assignment operator overload
                BaseEngine& operator= (BaseEngine&& base_engine);

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Public Member Variables *****/

                // TODO: ADD_FOR_BETA_RELEASE

            private:

                /***** Private Member Variables *****/

                // Local (i.e. Host) memory
                double* integral_buffer_;

                // Size of local buffer
                unsigned int buffer_size_;

                // Basis set pointer(s)
                std::vector<BasisSet*> basis_sets_;

                // Current integral operator
                IntegralOperator integral_operator_;

                // (Optional) Digestion functor for on-the-fly digestion of computed integrals
                IntegralConsumer integral_consumer_;

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Private Member Functions *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Virtual Private Member Functions *****/

                // Internal implementation of integral computation over a pair of shells
                virtual void compute_ints(unsigned int a_shell_index, unsigned int b_shell_index);

                // Internal implementation of integral computation over a quartet of shells
                virtual void compute_ints(unsigned int a_shell_index, unsigned int b_shell_index,
                                          unsigned int c_shell_index, unsigned int d_shell_index);

                // Internal implementation of integral computation over a triplet of shells (Reserved for RI-type computations)
                virtual void compute_ints(unsigned int a_shell_index, unsigned int b_shell_index, unsigned int c_shell_index);

                // Internal implementation of integral computation over a pair of shell sets
                virtual void compute_ints(unsigned int a_shell_set_index, unsigned int b_shell_set_index);

                // Internal implementation of integral computation over a triplet of shell sets (Reserved for RI-type computations)
                virtual void compute_ints(unsigned int a_shell_set_index, unsigned int b_shell_set_index, unsigned int c_shell_set_index);

                // Internal implementation of integral computation over a quartet of shell sets
                virtual void compute_ints(unsigned int a_shell_set_index, unsigned int b_shell_set_index,
                                          unsigned int c_shell_set_index, unsigned int d_shell_set_index);

                // Internal implementation of integral computation over a ShellPairSet
                virtual void compute_ints(ShellPairSet& shell_pair_set);

                // Internal implementation of integral computation over a ShellQuartetSet
                virtual void compute_ints(ShellQuartetSet& shell_quartet_set);

                // Internal implementation of integral computation over a ShellTupleSet
                virtual void compute_ints(ShellTupleSet& shell_tuple_set);

                // TODO: Add functionality for computing integrals over Shells/ShellSets/Shell(Pair/Quartet/Tuple)Sets in different basis sets

            };

        } // namespace libaccint::host::engine

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_ENGINE_BASE_ENGINE_HPP_ */