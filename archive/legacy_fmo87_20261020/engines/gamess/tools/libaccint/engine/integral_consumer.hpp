/*
 * integral_consumer.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_ENGINE_INTEGRAL_CONSUMER_HPP_
#define LIBACCINT_INCLUDE_ENGINE_INTEGRAL_CONSUMER_HPP_

// External header(s)

// Internal header(s)
#include "digestion_type.hpp"
#include "source/host/integral/digestion_params_base.hpp"

// STL header(s)

// Convenience typedef(s)

namespace libaccint {

    namespace host {

        namespace engine {

            /*! \Class IntegralConsumer
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief Function object for digesting integrals
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            class IntegralConsumer {

            public:

                /***** Constructor(s) and Destructor *****/

                // Default constructor
                IntegralConsumer();

                // Specialized constructor(s)
                // TODO

                // Deleted Copy constructor
                IntegralConsumer(const IntegralConsumer& int_consumer) = delete;

                // Deleted Move constructor
                IntegralConsumer(IntegralConsumer&& int_consumer) = delete;

                // Default destructor
                ~IntegralConsumer();

                /***** Public Member Functions *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Public Member Variables *****/

                // TODO: ADD_FOR_BETA_RELEASE

                /***** Operator Overloads *****/

                // Operator overload for particular digestion function
                void operator()(DigestionType digestion_type);

                // Deleted Copy assignment operator
                IntegralConsumer& operator= (const IntegralConsumer& int_consumer) = delete;

                // Deleted Move assignment operator
                IntegralConsumer& operator= (IntegralConsumer&& int_consumer) = delete;

                // TODO: Other operator overloads as needed

            private:

                /***** Private Member Functions *****/

                // TODO

                /***** Private Member Variables *****/

                // Pointer to base digestion parameters class
                // Specific parameters for a given digestion function are implemented in a 
                DigestionParamsBase* digestion_parameters_;

                // TODO: Etc.

            };

        } // namespace libaccint::host::engine

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_ENGINE_INTEGRAL_CONSUMER_HPP_ */
