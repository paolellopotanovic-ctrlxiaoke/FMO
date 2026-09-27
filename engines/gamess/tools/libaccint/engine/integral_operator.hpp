/*
 * integral_operator.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_ENGINE_INTEGRAL_OPERATOR_HPP_
#define LIBACCINT_INCLUDE_ENGINE_INTEGRAL_OPERATOR_HPP_

// External header(s)

// Internal header(s)
#include "integral_operator_type.hpp"
#include "source/host/integral/operator_params_base.hpp"

// STL header(s)

// Convenience typedef(s)

namespace libaccint {

    namespace host {

        namespace engine {

            /*! \Class IntegralOperator
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief Class for particular operator that integrals are evaluated over
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            class IntegralOperator {

            public:

                /***** Constructor(s) and Destructor *****/

                // Default constructor
                IntegralOperator();

                // Specialized constructor(s)
                // TODO

                // Deleted Copy constructor
                IntegralOperator(const IntegralOperator& integral_operator) = delete;

                // Deleted Move constructor
                IntegralOperator(IntegralOperator&& integral_operator) = delete;

                // Default destructor
                ~IntegralOperator();

                /***** Public Member Functions *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Public Member Variables *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Operator Overloads *****/

                // Deleted Copy assignment operator
                IntegralOperator& operator= (const IntegralOperator& integral_operator) = delete;

                // Deleted Move assignment operator
                IntegralOperator& operator= (IntegralOperator&& integral_operator) = delete;

                // TODO: Other operator overloads as needed

            private:

                /***** Private Member Functions *****/

                // TODO

                /***** Private Member Variables *****/

                // Specific operator type
                IntegralOperatorType int_operator_type_;

                // Pointer to base operator parameters class
                OperatorParamsBase* operator_parameters_;

                // TODO: Etc.

            };

        } // namespace libaccint::host::engine

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_ENGINE_INTEGRAL_OPERATOR_HPP_ */
