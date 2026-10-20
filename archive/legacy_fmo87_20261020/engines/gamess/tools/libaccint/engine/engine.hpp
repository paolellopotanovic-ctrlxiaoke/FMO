/*
 * engine.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_ENGINE_ENGINE_HPP_
#define LIBACCINT_INCLUDE_ENGINE_ENGINE_HPP_

// External header(s)

// Internal header(s)
#include "base_engine.hpp"
#include "compute_element.hpp"

// STL header(s)

namespace libaccint {

    namespace host {

        namespace engine {

            /*! \class Engine
            *  \ingroup TODO:ADD_GROUP
            *
            *  \brief Engine class template declaration
            *
            *  \author Bryce M. Westheimer
            *
            *  \date RELEASE_DATE
            *
            *  TODO: Etc.
            *
            */
            template<ComputeElement compute_element>
            class Engine : public BaseEngine;

        } // namespace libaccint::host::engine

    } // namespace libaccint::host

} // namespace libaccint

// Included template specialization headers, which themselves include implementation files
#include "source/host/engine/engine.hpp"    // Engine<CPU>
#include "source/device/engine/engine.hpp"  // Engine<GPU>

#endif /* LIBACCINT_INCLUDE_ENGINE_ENGINE_HPP_ */