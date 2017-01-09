/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#include "wfc.h"
#include "wfc_export.h"
#include <boost/python.hpp>
using namespace boost::python;


void export_wfc(){

  class_<wfc>("wfc",init<>())
    .def(init<wfc&,int,int,wfc&,int,int>())
    .def_readwrite("nspin",&wfc::nspin)
    .def_readwrite("gamma_only",&wfc::gamma_only)
    .def_readwrite("natoms",&wfc::natoms)
    .def_readwrite("tpiba",&wfc::tpiba)
    .def_readwrite("alat",&wfc::alat)
    .def_readwrite("omega",&wfc::omega)
    .def_readwrite("efermi",&wfc::efermi)
    .def_readwrite("cell_units",&wfc::cell_units)
    .def_readwrite("energy_units",&wfc::energy_units)
    .def_readwrite("nkpts",&wfc::nkpts)
    .def_readwrite("nbands",&wfc::nbands)
    .def_readwrite("npw",&wfc::npw)


    .def("QE_read_binary_wfc",&wfc::QE_read_binary_wfc)
    .def("QE_read_acsii_wfc",&wfc::QE_read_acsii_wfc)
    .def("QE_read_acsii_grid",&wfc::QE_read_acsii_grid)
    .def("QE_read_acsii_index",&wfc::QE_read_acsii_index)

    .def("complete",&wfc::complete)
    .def("normalize",&wfc::normalize)
    .def("restore",&wfc::restore)
    .def("set_latt_vectors",&wfc::set_latt_vectors)
    .def("set_reci_vectors",&wfc::set_reci_vectors)
    .def("compute_Hprime",&wfc::compute_Hprime)
  ;

  def("overlap",&overlap);
  def("nac",&nac);
  def("energy",&energy);
  def("ham",&ham);

}
