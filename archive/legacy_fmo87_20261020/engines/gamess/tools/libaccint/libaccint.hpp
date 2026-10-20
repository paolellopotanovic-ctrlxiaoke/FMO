/*
 * libaccint.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_HPP_
#define LIBACCINT_HPP_

// Core interface files
#include "core/atom.hpp"
#include "core/shell.hpp"
#include "core/shell_set.hpp"
#include "core/basis_set.hpp"

// Integral interface files
#include "integral/integral_operator.hpp"
#include "integral/shell_pair_set.hpp"
#include "integral/shell_quartet_set.hpp"
#include "integral/shell_tuple_set.hpp"

// Engine interface files
#include "engine/compute_element.hpp"
#include "engine/engine.hpp"

// Bootstrapping
#include "source/host/runtime/runtime.hpp"

namespace libaccint {

	// Using statements for necessary functions
	using libaccint::host::runtime::initialize;
	using libaccint::host::runtime::finalize;

	// Convenience typedef(s) and alias template(s)
	/***** Core *****/
	typedef libaccint::host::core::Atom		Atom;
	typedef libaccint::host::core::Shell	Shell;
	typedef libaccint::host::core::ShellSet	ShellSet;
	typedef libaccint::host::core::BasisSet BasisSet;
	
	/***** Integral *****/
	typedef libaccint::host::integral::IntegralOperator IntegralOperator;
	
	/***** Engine *****/
	typedef libaccint::host::engine::ComputeElement ComputeElement;
	template<ComputeEleement compute_element>
	using Engine = Engine<compute_element>;

	// For codes which intend to use each engine type separately
	typedef libaccint::host::engine::Engine<ComputeElement::CPU> HostEngine;
	typedef libaccint::host::engine::Engine<ComputeElement::GPU> DeviceEngine;

}
#endif /* LIBACCINT_HPP_ */