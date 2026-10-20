/*
 * Copyright (c) 2004-2005 The Trustees of Indiana University and Indiana
 *                         University Research and Technology
 *                         Corporation.  All rights reserved.
 * Copyright (c) 2004-2005 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2004-2005 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 * Copyright (c) 2004-2005 The Regents of the University of California.
 *                         All rights reserved.
 * Copyright (c) 2011 Cisco Systems, Inc.  All rights reserved.
 * Copyright (c) 2016      Research Organization for Information Science
 *                         and Technology (RIST). All rights reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 *
 * This file should be included by any file that needs full
 * version information for the OPAL project
 */

#ifndef OPAL_VERSIONS_H
#define OPAL_VERSIONS_H

#define OPAL_MAJOR_VERSION 4
#define OPAL_MINOR_VERSION 1
#define OPAL_RELEASE_VERSION 7
#define OPAL_GREEK_VERSION "rc1"
#define OPAL_WANT_REPO_REV @OPAL_WANT_REPO_REV@
#define OPAL_REPO_REV "v4.1.5-176-g6d9519e4c3"
#ifdef OPAL_VERSION
/* If we included version.h, we want the real version, not the
   stripped (no-r number) verstion */
#undef OPAL_VERSION
#endif
#define OPAL_VERSION "4.1.7rc1"
#define OPAL_CONFIGURE_CLI " \'--prefix=/build-result/hpcx-v2.22.1-gcc-inbox-ubuntu24.04-cuda12-x86_64/ompi\' \'--with-libevent=internal\' \'--enable-mpi1-compatibility\' \'--without-xpmem\' \'--with-cuda=/hpc/local/oss/cuda12.6.3/ubuntu24.04\' \'--with-slurm\' \'--with-platform=contrib/platform/mellanox/optimized\' \'--with-hcoll=/build-result/hpcx-v2.22.1-gcc-inbox-ubuntu24.04-cuda12-x86_64/hcoll\' \'--with-ucx=/build-result/hpcx-v2.22.1-gcc-inbox-ubuntu24.04-cuda12-x86_64/ucx\' \'--with-ucc=/build-result/hpcx-v2.22.1-gcc-inbox-ubuntu24.04-cuda12-x86_64/ucc\'"

#endif
