/*
 * digestion_type.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_ENGINE_DIGESTION_TYPE_HPP_
#define LIBACCINT_INCLUDE_ENGINE_DIGESTION_TYPE_HPP_

 // External header(s)

 // Internal header(s)

 // STL header(s)


namespace libaccint {

    namespace host {

        namespace engine {

            /*! \enum DigestionType
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief Digestion function type enumerator
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            enum class DigestionType {
                FockMatrix,
                CoulombMatrix,
                ExchangeMatrix,
                RIMatrix
                // TODO: Add all others as implemented
            };

        } // namespace libaccint::host::engine

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_ENGINE_DIGESTION_TYPE_HPP_ */
