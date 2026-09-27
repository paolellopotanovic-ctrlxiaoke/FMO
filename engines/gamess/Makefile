#BASH allows the pipeline to return a fail if any command fails
#Since the $(COMP) command below outputs "module.f" in warnings and errors
#replacing this string with the path to the real source file will make special
#features of VIM and EMACS work.
SHELL=bash -o pipefail

ifneq ($(wildcard install.info),)
       GMS_PATH := $(shell /bin/tcsh -fc 'source install.info && echo $$GMS_PATH')
else
ifneq ($(wildcard ../install.info),)
       GMS_PATH := $(shell /bin/tcsh -fc 'source ../install.info && echo $$GMS_PATH')
else
       $(error Please run 'config' first, to set up GAMESS compiling information)
endif
endif

#
# Import variables from install.info
#
GMS_VERSION := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_VERSION')
GMS_BUILD_PATH := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_BUILD_DIR')
GMS_MATHLIB := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_MATHLIB')
GMS_OPENMP := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_OPENMP')
GMS_MSUAUTO := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_MSUAUTO')
GMS_HPCCHEM := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_HPCCHEM')
GMS_LIBXC := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_LIBXC')
GMS_MDI := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_MDI')
GMS_FC := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_FORTRAN')
GMS_TARGET := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_TARGET')
GMS_VM2 := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_VM2')
TINKER := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$TINKER')
VB2000 := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$VB2000')
XMVB := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$XMVB')
NEO := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$NEO')
NBO := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$NBO')
RISM := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$RISM')
GMS_OPENMP_OFFLOAD := $(shell /bin/tcsh -fc 'source $(GMS_PATH)/install.info && echo $$GMS_OPENMP_OFFLOAD')
#
#   Choose C compiler for LibXC
#
GMS_CC := cc
ifeq ($(GMS_FC),gfortran)
        GMS_CC := gcc
endif
ifeq ($(GMS_FC),xlf)
        GMS_CC := xlc
endif
ifeq ($(GMS_TARGET),fj-a64fx)
ifeq ($(GMS_FC),frt)
        GMS_CC := fcc
else
ifeq ($(GMS_FC),frtpx)
        GMS_CC := fccpx
endif
endif
endif
ifeq ($(GMS_TARGET),ibm-bg)
        GMS_CC := bgxlc_r
endif
ifeq ($(GMS_TARGET),linux64)
ifeq ($(GMS_FC),pgfortran)
        GMS_CC := pgcc
else
ifeq ($(GMS_FC),ifort)
        GMS_CC := icc
else
        GMS_CC := gcc
endif
endif
endif
ifeq ($(GMS_TARGET),cray-xc)
        GMS_CC := cc
endif
ifeq ($(GMS_TARGET),aarch64)
        GMS_CC := gcc
endif
ifeq ($(GMS_TARGET),armflang)
        GMS_CC := armclang
endif

#
#   Adjust Fortran compiler name for CMake builds
#
ifeq ($(GMS_FC),oneapi-ifort)
        GMS_FC := ifort
endif
ifeq ($(GMS_FC),oneapi-ifx)
        GMS_FC := ifx
endif
#
#   LAPACK
#
USE_CMAKE := true
LAPACK_CMAKE_Fortran_FLAGS := ""
LAPACK_MAKE_INC := ""
ifeq ($(GMS_FC),xlf)
	LAPACK_CMAKE_Fortran_FLAGS := -DCMAKE_Fortran_FLAGS="-O3 -q64 -qintsize=8 -qthreaded -qrecur -qnosave -qstrict"
	LAPACK_MAKE_INC := make.inc.xlf.64
endif
ifeq ($(GMS_FC),gfortran)
	LAPACK_MAKE_INC := make.inc.gfortran.64
endif
ifeq ($(GMS_FC),nvfortran)
	LAPACK_MAKE_INC := make.inc.nvfortran.64
endif
ifeq ($(GMS_FC),ifort)
	LAPACK_MAKE_INC := make.inc.ifort.64
	ifeq ($(GMS_TARGET),win64)
		LAPACK_MAKE_INC := make.inc.ifort.win64
		USE_CMAKE := false
	endif
endif
ifeq ($(GMS_FC),ifx)
	LAPACK_MAKE_INC := make.inc.ifx.64
	USE_CMAKE := false
endif
ifeq ($(GMS_FC),aocc)
	LAPACK_MAKE_INC := make.inc.aocc.64
	GMS_FC := flang
	GMS_CC := clang
endif
#
undefine GMS_TARGET
#
GMS_3RD_PATH= $(GMS_PATH)/3rd-party
GMS_TOOLS_PATH= $(GMS_PATH)/tools

COMP       = $(GMS_PATH)/comp
COMPALL    = $(GMS_PATH)/compall
LKED       = $(GMS_PATH)/lked
OBJ_DIR    = $(GMS_BUILD_PATH)/object
OBJ        = $(patsubst %,$(OBJ_DIR)/%.o,$(SOURCEFILES))
LIS        = $(patsubst %,$(OBJ_DIR)/%.lis,$(SOURCEFILES))
ACT        = $(patsubst %,$(OBJ_DIR)/%.[fF],$(SOURCEFILES))
OBJ_MOD    = $(patsubst %,$(OBJ_DIR)/%.o,$(MODULES))
OBJ_MOD_SRC= $(patsubst %,$(OBJ_DIR)/%.F90,$(MODULES))
SRC_DIR    = $(GMS_PATH)/source
DDI_DIR    = $(GMS_BUILD_PATH)/ddi
DDI_LIB    = $(DDI_DIR)/libddi.a
DDIKICK    = $(GMS_BUILD_PATH)/ddikick.x
DDITEST    = $(GMS_BUILD_PATH)/ddi_test.x
DATASERVER = $(GMS_BUILD_PATH)/data_server.x
EXE        = gamess
GAMESS     = $(GMS_BUILD_PATH)/$(EXE).$(GMS_VERSION).x
CTAGS      = ctagsf
FTNCHEK    = ftnchek
FTNCHEK_OPT= -noextern -common=2 -arguments=1 -array=0 -calltree=0 -nocrossref -portability=tab -f77=automatic-array,relops,cycle-exit,case-construct -nonovice
FTNCHEK_OPT2= -noextern -common=2 -arguments=1 -array=0 -calltree=0 -nocrossref -resources=yes -f77=automatic-array,relops -portability=tab -nonovice

LIBXC_FLAGS  = -DDISABLE_VXC=OFF -DDISABLE_FXC=OFF -DDISABLE_KXC=OFF -DDISABLE_LXC=ON -DENABLE_FORTRAN=ON -DENABLE_XHOST=OFF -DBUILD_TESTING=ON
CMAKE        = $(shell if [[ `command -v cmake3` != "" ]]; then echo cmake3; else echo cmake; fi)

MDI_FLAGS  = -Dlibtype=STATIC -Dlanguage=Fortran -Dplugins=OFF

SEDSRC     = 's%\([a-z0-9]*\)\.f%$(SRC_DIR)\/\1\.src%'
SEDSRC2    = 's%$(OBJ_DIR)/\([a-z0-9]*\)\.f%$(SRC_DIR)/\1\.src%'

# 'parse_section' parses a file $(2) for
# "$(1)-start-section" and # "$(1)-end-section" markers.
# Between these markers it will look for lines containing "./comp".
# For each of these lines, the beginning of the line up until "./comp "
# will be removed and the result will be printed.
parse_section = $(shell sed -n '/$(1)-start-section/,/$(1)-end-section/{\:./comp:{s:^\s*./comp ::;p}}' $(2))

MODULES    = $(call parse_section,gamess-modules,$(COMPALL))

#All SOURCEFILES will be read from $(COMPALL)
#   this reads all lines between "source files" and top of TINKER plug-in
#   so what it picks out is the standard GAMESS list (no plug-ins included)
#   VB2000 is in a nonstandard place, so append its files manually,
#   it will be at the discretion of 'lked' to use the real or dummy code.
SOURCEFILES := $(call parse_section,gamess-base,$(COMPALL))

#   This is how we deal with optional builds
ifeq ($(GMS_MSUAUTO),true)
	MODULES := $(MODULES) $(call parse_section,msuauto-modules-enabled,$(COMPALL))
	SOURCEFILES := $(SOURCEFILES) $(call parse_section,msuauto-enabled,$(COMPALL))
else
	SOURCEFILES := $(SOURCEFILES) $(call parse_section,msuauto-disabled,$(COMPALL))
endif



ifeq ($(GMS_HPCCHEM),true)
	MODULES := $(MODULES) $(call parse_section,hpcchem-modules,$(COMPALL))
	SOURCEFILES := $(SOURCEFILES) $(call parse_section,hpcchem,$(COMPALL))
endif
ifeq ($(GMS_OPENMP),true)
	MODULES := $(MODULES) $(call parse_section,openmp-modules,$(COMPALL))
	SOURCEFILES := $(SOURCEFILES) $(call parse_section,openmp,$(COMPALL))
endif

ifeq ($(GMS_OPENMP_OFFLOAD),true)
	MODULES := $(MODULES) $(shell sed -n '/openmp-offload-modules-start-section/,/openmp-offload-modules-end-section/p' $(COMPALL) | grep "\.\/comp" | awk '{printf "%s ",$$2}')
	SOURCEFILES := $(SOURCEFILES) $(shell sed -n '/openmp-offload-start-section/,/openmp-offload-end-section/p' $(COMPALL) | grep "\.\/comp" | awk '{printf "%s ",$$2}')
endif

ifeq ($(GMS_LIBXC),true)
	MODULES := $(call parse_section,libxc-enabled,$(COMPALL)) $(MODULES)
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,libxc-enabled,$(COMPALL))
else
	MODULES := $(call parse_section,libxc-disabled,$(COMPALL)) $(MODULES)
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,libxc-disabled,$(COMPALL))
endif

ifeq ($(GMS_MDI),true)
	MODULES := $(MODULES) $(call parse_section,mdi-enabled,$(COMPALL))
else
	MODULES := $(MODULES) $(call parse_section,mdi-disabled,$(COMPALL))
endif

ifeq ($(GMS_VM2),false)
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,vm2-disabled,$(COMPALL))
endif

ifeq ($(TINKER),true)
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,tinker,$(COMPALL))
endif

ifeq ($(VB2000),true)
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,vb2000-enabled,$(COMPALL))
else
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,vb2000-disabled,$(COMPALL))
endif

ifeq ($(XMVB),true)
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,xmvb-enabled,$(COMPALL))
else
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,xmvb-disabled,$(COMPALL))
endif

ifeq ($(NEO),true)
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,neo-enabled,$(COMPALL))
else
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,neo-disabled,$(COMPALL))
endif

ifeq ($(NBO),false)
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,nbo-disabled,$(COMPALL))
endif

ifeq ($(RISM),true)
	MODULES := $(MODULES) $(call parse_section,rism-modules-enabled,$(COMPALL))
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,rism-enabled,$(COMPALL))
else
	MODULES := $(MODULES) $(call parse_section,rism-modules-disabled,$(COMPALL))
	SOURCEFILES_OPTIONAL := $(SOURCEFILES_OPTIONAL) $(call parse_section,rism-disabled,$(COMPALL))
endif
#$(info =============================)
#$(info SOURCEFILES is $(SOURCEFILES))
#$(info =============================)
#$(info SOURCEFILES_OPTIONAL is $(SOURCEFILES_OPTIONAL))
#$(info =============================)
#$(info MODULES is $(MODULES))
#$(info =============================)

SOURCEFILES := $(SOURCEFILES) $(SOURCEFILES_OPTIONAL)

#GAMESS executable

.PHONY : gamess
ifeq ($(GMS_MATHLIB),none)
gamess : | $(GMS_3RD_PATH)/lib/liblapack64.a
endif
ifeq ($(origin NPY),environment)
ifeq ($(NPY),true)
gamess : $(OBJ_DIR)/npy.o
endif
endif
gamess : | $(OBJ_DIR) $(GAMESS)
$(GAMESS) : $(DDI_LIB) $(OBJ_MOD) $(OBJ)
	@$(LKED) $(EXE) $(GMS_VERSION)

#Creat the object directory
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

#Modules that need to be compiled first
.PHONY : modules
modules : | $(OBJ_DIR) $(OBJ_MOD)
ifeq ($(origin NPY),environment)
ifeq ($(NPY),true)
modules : $(OBJ_DIR)/npy.o
endif
endif

$(OBJ_DIR)/neo.o :$(SRC_DIR)/../qmnuc/neo/neo.src
	$(call remove_single_file_silent,$(OBJ_DIR)/$(notdir $(basename $@)).f)
	@$(COMP) $(notdir $(basename $@)) 1>/dev/null # | sed -e $(SEDSRC)

$(OBJ_DIR)/neo%.o :$(SRC_DIR)/../qmnuc/neo/neo%.src
	$(call remove_single_file_silent,$(OBJ_DIR)/$(notdir $(basename $@)).f)
	@$(COMP) $(notdir $(basename $@)) 1>/dev/null # | sed -e $(SEDSRC)

$(OBJ_DIR)/ccsd3a%.o :$(SRC_DIR)/../msuauto/src/ccsd3a%.src
	$(call remove_single_file_silent,$(OBJ_DIR)/$(notdir $(basename $@)).f)
	@$(COMP) $(notdir $(basename $@)) 1>/dev/null # | sed -e $(SEDSRC)
sourcefiles : $(OBJ)
$(OBJ_DIR)/deaeom%.o :$(SRC_DIR)/../msuauto/src/deaeom%.src
	$(call remove_single_file_silent,$(OBJ_DIR)/$(notdir $(basename $@)).f)
	@$(COMP) $(notdir $(basename $@)) 1>/dev/null # | sed -e $(SEDSRC)
sourcefiles : $(OBJ)
$(OBJ_DIR)/eomdip%.o :$(SRC_DIR)/../msuauto/src/eomdip%.src
	$(call remove_single_file_silent,$(OBJ_DIR)/$(notdir $(basename $@)).f)
	@$(COMP) $(notdir $(basename $@)) 1>/dev/null # | sed -e $(SEDSRC)
sourcefiles : $(OBJ)

$(OBJ_DIR)/rsm_%.o :$(SRC_DIR)/../rsmcsed/src/rsm_%.src
	$(call remove_single_file_silent,$(OBJ_DIR)/$(notdir $(basename $@)).f)
	@$(COMP) $(notdir $(basename $@)) 1>/dev/null # | sed -e $(SEDSRC)

#Compile any module (I can't quite do "make $(module)" yet)
.PHONY : sourcefiles
sourcefiles : $(OBJ)

# Automatic prerequisite generation to enable builds with `make -j 4`
define OBJ_template =
$$(OBJ_DIR)/$(1).o: $$(OBJ_DIR)/$(2).o
endef
# Deduplicate modules
MODULES_DD:=$(shell r=""; for i in $(MODULES); do if ! [[ " $$r " == *" $$i "* ]]; then r="$$r $$i"; fi; done; echo $$r)
# Deduplicate sourcefiles, ignore sourcefiles that are already declared as modules
SOURCEFILES_DD:=$(shell r="";s="$(MODULES_DD)"; for i in $(SOURCEFILES); do if ! [[ " $$r " == *" $$i "* ]] && ! [[ " $$s " == *" $$i "* ]]; then r="$$r $$i"; fi; done; echo $$r)
# Each module depends on its predecessor; This forces the modules to be compiled one after another.
$(foreach idx,$(shell seq $$(($(words $(MODULES_DD)) - 1))),$(eval $(call OBJ_template,$(word $(shell echo $$(( $(idx)+1 )) ),$(MODULES_DD)),$(word $(idx),$(MODULES_DD)))))
# Each sourcefile depends on the last module; Together with the dependencies between modules
# this forces all modules to be compiled ahead of all sourcefiles
$(foreach sf,$(SOURCEFILES_DD),$(eval $(call OBJ_template,$(sf),$(lastword $(MODULES_DD)))))

$(OBJ_DIR)/%.o :$(SRC_DIR)/%.src
	$(call remove_single_file_silent,$(OBJ_DIR)/$(notdir $(basename $@)).f)
	@$(COMP) $(notdir $(basename $@)) 1>/dev/null # | sed -e $(SEDSRC)

#I like to keep the activated source in the OBJ_DIR for use with FTNCHEK
$(OBJ_DIR)/%.f :$(SRC_DIR)/%.src
	$(call remove_single_file_silent,$(OBJ_DIR)/$(notdir $(basename $@)).f)
	@$(COMP) $(notdir $(basename $@)) true
$(OBJ_DIR)/%.lis :$(SRC_DIR)/%.src $(OBJ_DIR)/%.f
	$(FTNCHEK) $(FTNCHEK_OPT) $(OBJ_DIR)/$(notdir $(basename $@)).f 2>&1 | sed -e $(SEDSRC2) > $@

# VB2000 is not detected as part of GAMESS, it is still set up as a plugin,
# located in non-standard places, so we need special rules to test their
# creation dates.
$(OBJ_DIR)/mod_vb2000.o : $(SRC_DIR)/../vb2000/SRC/mod_vb2000.src
	$(call remove_single_file_silent,$(OBJ_DIR)/mod_vb2000.f)
	@$(COMP) mod_vb2000 2>&1 1>/dev/null | sed -e $(SEDSRC)
$(OBJ_DIR)/vb2000.o : $(SRC_DIR)/../vb2000/SRC/vb2000.src
	$(call remove_single_file_silent,$(OBJ_DIR)/vb2000.f)
	@$(COMP) vb2000 2>&1 1>/dev/null | sed -e $(SEDSRC)
$(OBJ_DIR)/vb2gms.o : $(SRC_DIR)/../vb2000/SRC/vb2gms.src
	$(call remove_single_file_silent,$(OBJ_DIR)/vb2gms.f)
	@$(COMP) vb2gms 2>&1 1>/dev/null | sed -e $(SEDSRC)
$(OBJ_DIR)/vbdum.o : $(SRC_DIR)/../misc/vbdum.src
	$(call remove_single_file_silent,$(OBJ_DIR)/vbdum.f)
	@$(COMP) vbdum 2>&1 1>/dev/null | sed -e $(SEDSRC)
# and the .f for ftnchek
$(OBJ_DIR)/mod_vb2000.f : $(SRC_DIR)/../vb2000/SRC/mod_vb2000.src
	$(call remove_single_file_silent,$(OBJ_DIR)/mod_vb2000.f)
	@$(COMP) mod_vb2000 2>&1 1>/dev/null true | sed -e $(SEDSRC)
$(OBJ_DIR)/vb2000.f : $(SRC_DIR)/../vb2000/SRC/vb2000.src
	$(call remove_single_file_silent,$(OBJ_DIR)/vb2000.f)
	@$(COMP) vb2000 2>&1 1>/dev/null true | sed -e $(SEDSRC)
$(OBJ_DIR)/vb2gms.f : $(SRC_DIR)/../vb2000/SRC/vb2gms.src
	$(call remove_single_file_silent,$(OBJ_DIR)/vb2gms.f)
	@$(COMP) vb2gms 2>&1 1>/dev/null true | sed -e $(SEDSRC)
$(OBJ_DIR)/vbdum.f : $(SRC_DIR)/../misc/vbdum.src
	$(call remove_single_file_silent,$(OBJ_DIR)/vbdum.f)
	@$(COMP) vbdum 2>&1 1>/dev/null true | sed -e $(SEDSRC)
# and the .lis files for ftnchek
$(OBJ_DIR)/mod_vb2000.lis : $(OBJ_DIR)/mod_vb2000.f
	$(call remove_single_file_silent,$(OBJ_DIR)/mod_vb2000.f)
	@$(COMP) mod_vb2000 2>&1 1>/dev/null true | sed -e $(SEDSRC)
$(OBJ_DIR)/vb2000.lis : $(OBJ_DIR)/vb2000.f
	$(call remove_single_file_silent,$(OBJ_DIR)/vb2000.f)
	@$(COMP) vb2000 2>&1 1>/dev/null true | sed -e $(SEDSRC)
$(OBJ_DIR)/vb2gms.lis : $(OBJ_DIR)/vb2gms.f
	$(call remove_single_file_silent,$(OBJ_DIR)/vb2gms.f)
	@$(COMP) vb2gms 2>&1 1>/dev/null true | sed -e $(SEDSRC)
$(OBJ_DIR)/vbdum.lis : $(OBJ_DIR)/vbdum.f
	$(call remove_single_file_silent,$(OBJ_DIR)/vbdum.f)
	@$(COMP) vbdum 2>&1 1>/dev/null true | sed -e $(SEDSRC)

# NPY-for-Fortran-barebone
$(OBJ_DIR)/npy.o : | $(OBJ_DIR)
	cd $(OBJ_DIR) && gfortran -c -ffree-form -ffree-line-length-0 -std=f2003 -DNPY_UNIT=100 $(GMS_TOOLS_PATH)/NPY-for-Fortran-barebone/src/npy.F90

#I don't yet try to micro manage DDI compilation.
#Compilation is very quick anyway.
.PHONY : ddi
ddi : $(DDI_LIB) $(DDIKICK)
$(DDI_LIB) $(DDIKICK) :
	$(GMS_PATH)/ddi/compddi
	@if [ -e $(GMS_PATH)/ddi/ddikick.x ]; then mv $(GMS_PATH)/ddi/ddikick.x $(DDIKICK);fi
	@if [ -e $(GMS_PATH)/ddi/ddi_test.x ]; then mv $(GMS_PATH)/ddi/ddi_test.x $(DDITEST);fi
	@if [ -e $(GMS_PATH)/ddi/data_server.x ]; then mv $(GMS_PATH)/ddi/data_server.x $(DATASERVER);fi
	@if [ -e $(GMS_PATH)/ddi/libddi.a ] && ! [ "$(GMS_PATH)/ddi/libddi.a" -ef "$(DDI_LIB)" ]; then mv $(GMS_PATH)/ddi/libddi.a  $(DDI_LIB);fi


#### Auxillary targets ####

#at runtime, extra flags can be added (like -e to make emacs compitable tags)
.PHONY : tags
tags : $(SRC_DIR)/*.src
	$(CTAGS) $(CTAGS_FLAGS) --langmap=fortran:.src $(SRC_DIR)/* $(GMS_PATH)/ddi/src/*

#FTNCHEK: a summary is printed to STDOUT along with the gamess.lis and individual lis files in OBJ_DIR
gamess.lis : $(ACT) $(LIS)
	$(FTNCHEK) $(FTNCHEK_OPT2) $(ACT) 2>&1 | sed -e $(SEDSRC2) > gamess.lis
.PHONY : check
check : gamess.lis
	@awk '{if($$0 ~ "[1-9]+ syntax errors") print "------\n",$$0,"\n------"}' gamess.lis
	@awk '{if($$0 ~ "[1-9]+ warning[s]* issued") print $$0}' gamess.lis

#Build a tarball based on a git commit or tag
COMMIT=HEAD
TAR_VERSION=$(GMS_VERSION)

.PHONY : dist
dist :
	@echo "Creating $(GMS_BUILD_PATH)/gamess-$(TAR_VERSION).tar.gz for commit $(COMMIT) in branch `git name-rev --name-only $(COMMIT)`"
	@cd $(GMS_PATH) && git archive --format=tar --prefix=gamess-$(TAR_VERSION)/ $(COMMIT) | gzip >$(GMS_BUILD_PATH)/gamess-$(TAR_VERSION).tar.gz

.PHONY : clean
clean :
	$(call remove_multiple_files,$(OBJ))
	$(call remove_multiple_files,$(OBJ_DIR)/*.mod)
	$(call remove_multiple_files,$(OBJ_MOD))
	$(call remove_multiple_files,$(OBJ_MOD_SRC))
	$(call remove_multiple_files,$(ACT))
	$(call remove_multiple_files,$(LIS))
	$(call remove_single_file,$(GAMESS))
	$(call remove_single_file,gamess.lis)

.PHONY : clean_ddi
clean_ddi :
	$(call remove_single_file,$(DDI_LIB))
	$(call remove_single_file,$(DDIKICK))
	$(call remove_single_file,$(DDITEST))
	$(call remove_single_file,$(DATASERVER))

-PHONY : cleaned_gamess
cleaned_gamess: clean clean_ddi
	+@$(MAKE)

############ custom file removal functions #####
# In order to avoid 'rm -f' calls the following functions were created.
# Before deleting a path they check, whether it is really a file and
# whether the user has write permissions on this file. This avoids
# error messages or prompts from 'rm'.
#
# 'remove_single_file' has a single argument - the path of the file
# that should be deleted. It prints this name out to stdout.
#
# 'remove_single_file_silent' has a single argument - the path of the
# file that should be deleted. It is silent.
#
# 'remove_multiple_files' has a single argument - the list of file paths
# that should be deleted. It prints to stdout and uses one a separate
# line for each path.
#
#

define remove_single_file
    @echo "rm $(1)" ; \
    if [[ -f "$(1)" && -w "$(1)" ]]; then \
        rm "$(1)" ; \
    fi
endef

define remove_single_file_silent
    @if [[ -f "$(1)" && -w "$(1)" ]]; then \
        rm "$(1)" ; \
    fi
endef

define remove_multiple_files
    @for file in $(1) ;  do \
        echo "rm $$file"; \
        if [[ -f "$$file" && -w "$$file" ]]; then \
            rm "$$file"; \
        fi; \
    done
endef

############ The new CHECKTEST ##################
# Exam inputs are divided into subdirectory collections.  The inputs in the
# "standard" directory are default.
# More directories can be added to easily run all tests of a certain classification
# such as "mp2" or "gradient".  Soft links can be used if one wishes to place the same
# input file in two classifications.
# To run a different subdirectory, specify the name on the command line
#
# make checktest EXAM_PATH=mp2
#
# The target "checktest" will run the exam inputs and perform the error checking.
# If an examXX.log file exists, it will not be rerun but only check.  If the
# user wishes to only a few of the exams with in the chosen $(EXAM_PATH):
#
# make checktest EXAMS=exam05
# make checktest EXAMS="exam05 exam30"
#
# To add new exam, two files are required.  The chkabs script in the standard
# subdirectory is used by default.  It is recommended that names of your input
# files do not exceed 10 characters.  Any further characters will be dropped from
# the output of make checktest.
#
# The input
#    $(GMS_PATH)/tests/$(EXAM_PATH)/mytest.inp
# The script (use already standard/scripts/examXX scripts as templates)
#    $(GMS_PATH)/tests/$(EXAM_PATH)/scripts/mytest
#
# The new exam can be run and tested by itself
#    make checktest EXAM_PATH=mytests EXAMS=mytest
#
#

RUNGMS=$(GMS_PATH)/tools/localgms
NCPUS=1

EXAM_PATH = standard

GMS_EXT = log

EXAM_INPUTS := $(wildcard $(GMS_PATH)/tests/$(EXAM_PATH)/*.inp)
EXAMS = $(patsubst $(GMS_PATH)/tests/$(EXAM_PATH)/%.inp,%,$(EXAM_INPUTS))

EXAM_OUTPUT_DIR = $(GMS_BUILD_PATH)/tests/$(EXAM_PATH)
EXAM_LOG = $(patsubst %,$(EXAM_OUTPUT_DIR)/%.$(GMS_EXT),$(EXAMS))
EXAM_RST = $(patsubst %,$(EXAM_OUTPUT_DIR)/%.rst,$(EXAMS))
EXAM_DAT = $(patsubst %,$(EXAM_OUTPUT_DIR)/%.dat,$(EXAMS))
EXAM_TRJ = $(patsubst %,$(EXAM_OUTPUT_DIR)/%.trj,$(EXAMS))
EXAM_CIM = $(patsubst %,$(EXAM_OUTPUT_DIR)/%.cim,$(EXAMS))
EXAM_DMN = $(patsubst %,$(EXAM_OUTPUT_DIR)/%.dmn,$(EXAMS))
#This is a dummy file that is never created, but allows all tests to be checked all the time.
EXAM_CHK = $(patsubst %,$(EXAM_OUTPUT_DIR)/%.chk,$(EXAMS))

export GMSPATH=$(GMS_PATH)
export EXEPATH=$(GMS_BUILD_PATH)

$(EXAM_OUTPUT_DIR)/%.$(GMS_EXT) : $(GAMESS) $(GMS_PATH)/tests/$(EXAM_PATH)/%.inp
	@mkdir -p $(EXAM_OUTPUT_DIR)
	@cp $(GMS_PATH)/tests/$(EXAM_PATH)/$(notdir $(basename $@)).inp $(EXAM_OUTPUT_DIR)/$(notdir $(basename $@)).F05
	$(call remove_single_file_silent,$(EXAM_OUTPUT_DIR)/$(notdir $(basename $@)).dat)
	$(call remove_single_file_silent,$(EXAM_OUTPUT_DIR)/$(notdir $(basename $@)).rst)
	$(call remove_single_file_silent,$(EXAM_OUTPUT_DIR)/$(notdir $(basename $@)).trj)
	$(call remove_single_file_silent,$(EXAM_OUTPUT_DIR)/$(notdir $(basename $@)).cim)
	$(call remove_single_file_silent,$(EXAM_OUTPUT_DIR)/$(notdir $(basename $@)).dmn)
	$(call remove_single_file_silent,$(EXAM_OUTPUT_DIR)/$(notdir $(basename $@)).$(GMS_EXT))
	cd $(EXAM_OUTPUT_DIR) && $(RUNGMS) $(notdir $(basename $@)) $(GMS_VERSION) $(NCPUS) > $(notdir $(basename $@)).$(GMS_EXT) 2>&1

.PHONY : exams
exams : $(EXAM_LOG)


#These environment variables are required to run the
# the checking scripts in $(GMS_PATH)/tests/$(EXAM_PATH)/scripts
#Tolerances can be changed on the command line
#  make checktest tolE=6e-5
export tolC=0.3
export tolD=0.0001
export tolE=0.00000001
export tolG=0.00001
export tolH=0.0001
export tolI=0.0001
export tolL=0.1
export tolO=0.0001
export tolP=0.0001
export tolR=0.0001
export tolS=0.01
export tolT=0.000001
export tolV=0.00000001
export tolW=0.1
export tolX=0.00001
export tolIND=0.000001
export tolFREQ=0.4
#If using chkabs, only one needs to be supplied for all tests
export chkabs=$(GMS_PATH)/tests/standard/scripts/chkabs

CHK_PATH=$(GMS_PATH)/tests/$(EXAM_PATH)/scripts

#print the restuls of the examXX script
CHK_AWK = -v exam=$(notdir $(basename $@)) '{ printf "%-10s %-63s",exam,$$0 }';if [ $$? != 0 ]; then echo "!!Failed";else echo Passed;fi
#determine that the exam finished
GREP_TERM = grep -q "TERMINATED NORMALLY" $(EXAM_OUTPUT_DIR)/$(notdir $(basename $@)).$(GMS_EXT);if [ $$? != 0 ]; then echo "$(notdir $(basename $@))     !!! TERMINATED ABNORMALLY !!!";exit;fi

#check in single exam
$(EXAM_OUTPUT_DIR)/%.chk : $(EXAM_OUTPUT_DIR)/%.$(GMS_EXT)
	@cd $(EXAM_OUTPUT_DIR) && $(GREP_TERM) && $(CHK_PATH)/$(notdir $(basename $@)) $(notdir $(basename $@)).$(GMS_EXT) | awk $(CHK_AWK)

#since .chk files are never created make will always check all entries in EXAMS.
.PHONY : test checktest
test checktest: $(EXAM_LOG) $(EXAM_CHK)

.PHONY : clean_exams
clean_exams :
	$(call remove_multiple_files,$(EXAM_RST))
	$(call remove_multiple_files,$(EXAM_DAT))
	$(call remove_multiple_files,$(EXAM_LOG))
	$(call remove_multiple_files,$(EXAM_TRJ))
	$(call remove_multiple_files,$(EXAM_CIM))
	$(call remove_multiple_files,$(EXAM_DMN))

.PHONY: libxc libmdi

ifeq ($(GMS_LIBXC),true)
gamess : | $(GMS_3RD_PATH)/lib/libxc.a
endif

libxc: $(GMS_3RD_PATH)/lib/libxc.a

$(OBJ_DIR)/functionals.o : $(OBJ_DIR)/xc_f03_lib_m.mod

$(OBJ_DIR)/xc_f03_lib_m.mod : $(GMS_3RD_PATH)/include/xc_f03_lib_m.mod
	cp $(GMS_3RD_PATH)/include/xc_f03_lib_m.mod $(OBJ_DIR)

$(GMS_3RD_PATH)/include/xc_f03_lib_m.mod : $(GMS_3RD_PATH)/lib/libxc.a

$(GMS_3RD_PATH)/lib/libxc.a :
	cd $(GMS_3RD_PATH)/libxc && $(CMAKE) -H. -Bobjdir -DCMAKE_C_COMPILER=$(GMS_CC) -DCMAKE_Fortran_COMPILER=$(GMS_FC) -DCMAKE_INSTALL_PREFIX=$(GMS_3RD_PATH) $(LIBXC_FLAGS) --debug-output && cd objdir && $(MAKE) VERBOSE=1 && $(MAKE) test && $(MAKE) install
	if [[ ! -d $(GMS_3RD_PATH)/lib ]]; then ln -s $(GMS_3RD_PATH)/lib64 $(GMS_3RD_PATH)/lib; fi

libmdi : | $(OBJ_DIR)
	cd $(GMS_3RD_PATH)/mdi && $(CMAKE) -H. -Bobjdir -DCMAKE_C_COMPILER=$(GMS_CC) -DCMAKE_Fortran_COMPILER=$(GMS_FC) -DCMAKE_INSTALL_PREFIX=$(MDI_INSTALL) $(MDI_FLAGS) && cd objdir && $(MAKE) install
	cp $(GMS_3RD_PATH)/include/mdi/mdi.mod $(OBJ_DIR)
	if [[ ! -d $(GMS_3RD_PATH)/lib ]]; then ln -s $(GMS_3RD_PATH)/lib64 $(GMS_3RD_PATH)/lib; fi

ifeq ($(USE_CMAKE),true)
.PHONY : lapack
lapack : $(GMS_3RD_PATH)/lib/liblapack64.a

$(GMS_3RD_PATH)/lib/liblapack64.a :
	cd $(GMS_3RD_PATH)/lapack && $(CMAKE) -H. -Bobjdir -DCMAKE_C_COMPILER=$(GMS_CC) -DCMAKE_Fortran_COMPILER=$(GMS_FC) -DCMAKE_INSTALL_PREFIX=$(GMS_3RD_PATH) -DBUILD_INDEX64=ON -DBUILD_SINGLE=OFF -DBUILD_COMPLEX=OFF -DCBLAS=OFF -DLAPACKE=OFF $(LAPACK_CMAKE_Fortran_FLAGS) && cd objdir && $(MAKE) VERBOSE=1 && $(MAKE) install
	if [[ ! -d $(GMS_3RD_PATH)/lib ]]; then ln -s $(GMS_3RD_PATH)/lib64 $(GMS_3RD_PATH)/lib; fi
endif

ifeq ($(USE_CMAKE),false)
.PHONY : lapack
lapack : $(GMS_3RD_PATH)/lib/liblapack64.a

LAPACK_VERSION := 3.10.1

$(GMS_3RD_PATH)/lapack :
	mkdir -p $(GMS_3RD_PATH)/lapack

$(GMS_PATH)/lapack-$(LAPACK_VERSION).tar.gz : $(GMS_3RD_PATH)/lapack
	if [[ ! -f $(GMS_PATH)/lapack-$(LAPACK_VERSION).tar.gz ]]; then curl -LJO https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v$(LAPACK_VERSION).tar.gz ${SSL_REVOKE} -o $(GMS_PATH)/lapack-$(LAPACK_VERSION).tar.gz; fi
	tar -xzf $(GMS_PATH)/lapack-$(LAPACK_VERSION).tar.gz -C $(GMS_3RD_PATH)/lapack --strip-components=1
	patch -d ${GMS_PATH} -p1 < ${GMS_TOOLS_PATH}/lapack/lapacke.CMakeList.patch
	patch -d ${GMS_PATH} -p1 < ${GMS_TOOLS_PATH}/lapack/CheckLAPACKCompilerFlags.patch
	patch -d ${GMS_PATH} -p1 < ${GMS_TOOLS_PATH}/lapack/no_null_strings.source.patch

$(GMS_3RD_PATH)/lapack/make.inc : $(GMS_3RD_PATH)/lapack $(GMS_PATH)/lapack-$(LAPACK_VERSION).tar.gz
	cp $(GMS_TOOLS_PATH)/lapack/$(LAPACK_MAKE_INC) $(GMS_3RD_PATH)/lapack/make.inc

$(GMS_3RD_PATH)/lapack/librefblas.a : $(GMS_3RD_PATH)/lapack/make.inc $(GMS_PATH)/lapack-$(LAPACK_VERSION).tar.gz
	cd $(GMS_3RD_PATH)/lapack && $(MAKE) blaslib

$(GMS_3RD_PATH)/lapack/liblapack.a : $(GMS_3RD_PATH)/lapack/make.inc $(GMS_PATH)/lapack-$(LAPACK_VERSION).tar.gz
	cd $(GMS_3RD_PATH)/lapack && $(MAKE) lapacklib

$(GMS_3RD_PATH)/lib :
	mkdir -p $(GMS_3RD_PATH)/lib

ifeq ($(GMS_FC),xlf)
$(GMS_3RD_PATH)/lib/liblapack64.a : $(GMS_3RD_PATH)/lapack/liblapack.a $(GMS_3RD_PATH)/lib
else
$(GMS_3RD_PATH)/lib/liblapack64.a : $(GMS_3RD_PATH)/lapack/librefblas.a $(GMS_3RD_PATH)/lapack/liblapack.a $(GMS_3RD_PATH)/lib
endif
	if [[ -f $(GMS_3RD_PATH)/lapack/librefblas.a ]]; then ln -fs $(GMS_3RD_PATH)/lapack/librefblas.a $(GMS_3RD_PATH)/lib/libblas64.a; fi
	if [[ -f $(GMS_3RD_PATH)/lapack/liblapack.a ]]; then ln -fs $(GMS_3RD_PATH)/lapack/liblapack.a $(GMS_3RD_PATH)/lib/liblapack64.a; fi

.PHONY : clean_lapack
clean_lapack :
	rm -rf $(GMS_3RD_PATH)/lapack
	$(call remove_single_file,lapack-$(LAPACK_VERSION).tar.gz)
	$(call remove_single_file,$(GMS_3RD_PATH)/lib/libblas64.a)
	$(call remove_single_file,$(GMS_3RD_PATH)/lib/liblapack64.a)
endif
