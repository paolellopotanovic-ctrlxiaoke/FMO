/*
 * compute_element.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_ENGINE_COMPUTE_ELEMENT_HPP_
#define LIBACCINT_INCLUDE_ENGINE_COMPUTE_ELEMENT_HPP_

// External header(s)

// Internal header(s)

// STL header(s)


namespace libaccint {

    namespace host {

        namespace engine {

            /*! \enum ComputeElement
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief Compute element enumerator
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            enum class ComputeElement {
                CPU,
                GPU
                // TODO: Add all others as implemented
            };

        } // namespace libaccint::host::engine

    } // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_ENGINE_COMPUTE_ELEMENT_HPP_ */
