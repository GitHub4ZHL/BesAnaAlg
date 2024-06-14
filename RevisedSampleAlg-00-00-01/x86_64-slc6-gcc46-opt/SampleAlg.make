#-- start of make_header -----------------

#====================================
#  Library SampleAlg
#
#   Generated Wed Dec  9 13:13:50 2020  by leo591653959
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_SampleAlg_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_SampleAlg_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_SampleAlg

SampleAlg_tag = $(tag)

#cmt_local_tagfile_SampleAlg = $(SampleAlg_tag)_SampleAlg.make
cmt_local_tagfile_SampleAlg = $(bin)$(SampleAlg_tag)_SampleAlg.make

else

tags      = $(tag),$(CMTEXTRATAGS)

SampleAlg_tag = $(tag)

#cmt_local_tagfile_SampleAlg = $(SampleAlg_tag).make
cmt_local_tagfile_SampleAlg = $(bin)$(SampleAlg_tag).make

endif

include $(cmt_local_tagfile_SampleAlg)
#-include $(cmt_local_tagfile_SampleAlg)

ifdef cmt_SampleAlg_has_target_tag

cmt_final_setup_SampleAlg = $(bin)setup_SampleAlg.make
cmt_dependencies_in_SampleAlg = $(bin)dependencies_SampleAlg.in
#cmt_final_setup_SampleAlg = $(bin)SampleAlg_SampleAlgsetup.make
cmt_local_SampleAlg_makefile = $(bin)SampleAlg.make

else

cmt_final_setup_SampleAlg = $(bin)setup.make
cmt_dependencies_in_SampleAlg = $(bin)dependencies.in
#cmt_final_setup_SampleAlg = $(bin)SampleAlgsetup.make
cmt_local_SampleAlg_makefile = $(bin)SampleAlg.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)SampleAlgsetup.make

#SampleAlg :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'SampleAlg'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = SampleAlg/
#SampleAlg::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

SampleAlglibname   = $(bin)$(library_prefix)SampleAlg$(library_suffix)
SampleAlglib       = $(SampleAlglibname).a
SampleAlgstamp     = $(bin)SampleAlg.stamp
SampleAlgshstamp   = $(bin)SampleAlg.shstamp

SampleAlg :: dirs  SampleAlgLIB
	$(echo) "SampleAlg ok"

#-- end of libary_header ----------------

SampleAlgLIB :: $(SampleAlglib) $(SampleAlgshstamp)
	@/bin/echo "------> SampleAlg : library ok"

$(SampleAlglib) :: $(bin)JpsiToPhiEtaAlg.o $(bin)alg_load.o $(bin)alg_entries.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(SampleAlglib) $?
	$(lib_silent) $(ranlib) $(SampleAlglib)
	$(lib_silent) cat /dev/null >$(SampleAlgstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(SampleAlglibname).$(shlibsuffix) :: $(SampleAlglib) $(SampleAlgstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" SampleAlg $(SampleAlg_shlibflags)

$(SampleAlgshstamp) :: $(SampleAlglibname).$(shlibsuffix)
	@if test -f $(SampleAlglibname).$(shlibsuffix) ; then cat /dev/null >$(SampleAlgshstamp) ; fi

SampleAlgclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)JpsiToPhiEtaAlg.o $(bin)alg_load.o $(bin)alg_entries.o

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

ifeq ($(INSTALLAREA),)
installarea = $(CMTINSTALLAREA)
else
ifeq ($(findstring `,$(INSTALLAREA)),`)
installarea = $(shell $(subst `,, $(INSTALLAREA)))
else
installarea = $(INSTALLAREA)
endif
endif

install_dir = ${installarea}/${CMTCONFIG}/lib
SampleAlginstallname = $(library_prefix)SampleAlg$(library_suffix).$(shlibsuffix)

SampleAlg :: SampleAlginstall

install :: SampleAlginstall

SampleAlginstall :: $(install_dir)/$(SampleAlginstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(SampleAlginstallname) :: $(bin)$(SampleAlginstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(SampleAlginstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(SampleAlginstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(SampleAlginstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(SampleAlginstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(SampleAlginstallname) $(install_dir)/$(SampleAlginstallname); \
	      echo `pwd`/$(SampleAlginstallname) >$(install_dir)/$(SampleAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(SampleAlginstallname), no installation directory specified"; \
	  fi; \
	fi

SampleAlgclean :: SampleAlguninstall

uninstall :: SampleAlguninstall

SampleAlguninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(SampleAlginstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(SampleAlginstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(SampleAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(SampleAlginstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),SampleAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)JpsiToPhiEtaAlg.d

$(bin)$(binobj)JpsiToPhiEtaAlg.d :

$(bin)$(binobj)JpsiToPhiEtaAlg.o : $(cmt_final_setup_SampleAlg)

$(bin)$(binobj)JpsiToPhiEtaAlg.o : $(src)JpsiToPhiEtaAlg.cxx
	$(cpp_echo) $(src)JpsiToPhiEtaAlg.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(SampleAlg_pp_cppflags) $(lib_SampleAlg_pp_cppflags) $(JpsiToPhiEtaAlg_pp_cppflags) $(use_cppflags) $(SampleAlg_cppflags) $(lib_SampleAlg_cppflags) $(JpsiToPhiEtaAlg_cppflags) $(JpsiToPhiEtaAlg_cxx_cppflags)  $(src)JpsiToPhiEtaAlg.cxx
endif
endif

else
$(bin)SampleAlg_dependencies.make : $(JpsiToPhiEtaAlg_cxx_dependencies)

$(bin)SampleAlg_dependencies.make : $(src)JpsiToPhiEtaAlg.cxx

$(bin)$(binobj)JpsiToPhiEtaAlg.o : $(JpsiToPhiEtaAlg_cxx_dependencies)
	$(cpp_echo) $(src)JpsiToPhiEtaAlg.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(SampleAlg_pp_cppflags) $(lib_SampleAlg_pp_cppflags) $(JpsiToPhiEtaAlg_pp_cppflags) $(use_cppflags) $(SampleAlg_cppflags) $(lib_SampleAlg_cppflags) $(JpsiToPhiEtaAlg_cppflags) $(JpsiToPhiEtaAlg_cxx_cppflags)  $(src)JpsiToPhiEtaAlg.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),SampleAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)alg_load.d

$(bin)$(binobj)alg_load.d :

$(bin)$(binobj)alg_load.o : $(cmt_final_setup_SampleAlg)

$(bin)$(binobj)alg_load.o : $(src)components/alg_load.cxx
	$(cpp_echo) $(src)components/alg_load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(SampleAlg_pp_cppflags) $(lib_SampleAlg_pp_cppflags) $(alg_load_pp_cppflags) $(use_cppflags) $(SampleAlg_cppflags) $(lib_SampleAlg_cppflags) $(alg_load_cppflags) $(alg_load_cxx_cppflags) -I../src/components $(src)components/alg_load.cxx
endif
endif

else
$(bin)SampleAlg_dependencies.make : $(alg_load_cxx_dependencies)

$(bin)SampleAlg_dependencies.make : $(src)components/alg_load.cxx

$(bin)$(binobj)alg_load.o : $(alg_load_cxx_dependencies)
	$(cpp_echo) $(src)components/alg_load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(SampleAlg_pp_cppflags) $(lib_SampleAlg_pp_cppflags) $(alg_load_pp_cppflags) $(use_cppflags) $(SampleAlg_cppflags) $(lib_SampleAlg_cppflags) $(alg_load_cppflags) $(alg_load_cxx_cppflags) -I../src/components $(src)components/alg_load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),SampleAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)alg_entries.d

$(bin)$(binobj)alg_entries.d :

$(bin)$(binobj)alg_entries.o : $(cmt_final_setup_SampleAlg)

$(bin)$(binobj)alg_entries.o : $(src)components/alg_entries.cxx
	$(cpp_echo) $(src)components/alg_entries.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(SampleAlg_pp_cppflags) $(lib_SampleAlg_pp_cppflags) $(alg_entries_pp_cppflags) $(use_cppflags) $(SampleAlg_cppflags) $(lib_SampleAlg_cppflags) $(alg_entries_cppflags) $(alg_entries_cxx_cppflags) -I../src/components $(src)components/alg_entries.cxx
endif
endif

else
$(bin)SampleAlg_dependencies.make : $(alg_entries_cxx_dependencies)

$(bin)SampleAlg_dependencies.make : $(src)components/alg_entries.cxx

$(bin)$(binobj)alg_entries.o : $(alg_entries_cxx_dependencies)
	$(cpp_echo) $(src)components/alg_entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(SampleAlg_pp_cppflags) $(lib_SampleAlg_pp_cppflags) $(alg_entries_pp_cppflags) $(use_cppflags) $(SampleAlg_cppflags) $(lib_SampleAlg_cppflags) $(alg_entries_cppflags) $(alg_entries_cxx_cppflags) -I../src/components $(src)components/alg_entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: SampleAlgclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(SampleAlg.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

SampleAlgclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library SampleAlg
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)SampleAlg$(library_suffix).a $(library_prefix)SampleAlg$(library_suffix).s? SampleAlg.stamp SampleAlg.shstamp
#-- end of cleanup_library ---------------
