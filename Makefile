# DNSdata/Makefile

NAME := DNSdata
OBJS := sgrid_$(NAME).o DNSdata.o DNS_CTS.o setADMvars.o TOVstar.o \
	DNSgrid.o DNS_compute_new_q_instar.o DNS_compute_new_q_atXYZ.o \
 	DNS_set_restmassintegrand.o DNS_set_J_ADM_VolInt_integrand.o \
	DNS_set_M_ADM_VolInt_integrand.o DNS_set_P_ADM_VolInt_integrand.o \
	DNS_Interpolate_ADMvars.o \
	set_DNSdata_Sigma_BCs.o DNS_BCs.o DNS_set_dlnIntegEuler.o \
	DNS_compute_chi.o DNS_EoS.o

include $(TOP)/Makefile.subdirs

DNS_CTS.c: DNS_CTS.m
	math < DNS_CTS.m > /dev/null

setADMvars.c: setADMvars.m
	math < setADMvars.m > /dev/null

DNS_compute_new_q_instar.c: DNS_compute_new_q_instar.m
	math < DNS_compute_new_q_instar.m > /dev/null

DNS_compute_new_q_atXYZ.c: DNS_compute_new_q_atXYZ.m
	math < DNS_compute_new_q_atXYZ.m > /dev/null

DNS_set_restmassintegrand.c: DNS_set_restmassintegrand.m
	math < DNS_set_restmassintegrand.m > /dev/null

DNS_set_J_ADM_VolInt_integrand.c: DNS_set_J_ADM_VolInt_integrand.m
	math < DNS_set_J_ADM_VolInt_integrand.m > /dev/null

DNS_set_M_ADM_VolInt_integrand.c: DNS_set_M_ADM_VolInt_integrand.m
	math < DNS_set_M_ADM_VolInt_integrand.m > /dev/null

DNS_set_P_ADM_VolInt_integrand.c: DNS_set_P_ADM_VolInt_integrand.m
	math < DNS_set_P_ADM_VolInt_integrand.m > /dev/null

set_DNSdata_Sigma_BCs.c: set_DNSdata_Sigma_BCs.m
	math < set_DNSdata_Sigma_BCs.m > /dev/null

DNS_set_dlnIntegEuler.c: DNS_set_dlnIntegEuler.m
	math < DNS_set_dlnIntegEuler.m > /dev/null
