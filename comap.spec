%define _basename comap
%define _version 1.5.2
%define _release 1
%define _prefix /usr

URL: http://home.gna.org/comap/

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: Julien Dutheil
Source: http://biopp.univ-montp2.fr/repos/sources/comap/%{_basename}-%{_version}.tar.gz
Summary: The CoMap package
Group: Productivity/Scientific/Other

Requires: libbpp-phyl9 = 2.2.0
Requires: libbpp-seq9 = 2.2.0
Requires: libbpp-core2 = 2.2.0

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.6.0
BuildRequires: gcc-c++ >= 4.0.0
BuildRequires: groff
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core2 = 2.2.0
BuildRequires: libbpp-core-devel = 2.2.0
BuildRequires: libbpp-seq9 = 2.2.0
BuildRequires: libbpp-seq-devel = 2.2.0
BuildRequires: libbpp-phyl9 = 2.2.0
BuildRequires: libbpp-phyl-devel = 2.2.0


AutoReq: yes
AutoProv: yes

%if 0%{?mandriva_version}
%if %{mandriva_version} >= 2011
BuildRequires: xz
%define zipext xz
%else
BuildRequires: lzma
%define zipext lzma
%endif
%else
%if 0%{?distribution:1} && "%{distribution}" == "Mageia"
BuildRequires: xz
%define zipext xz
%else
#For all other distributions:
BuildRequires: gzip
%define zipext gz
%endif
%endif

%description
Includes programs:
 - CoMap, (co)substitution mapping and coevolution detection,
 - MICA, Mutual Information Coevolution Analysis.

%prep
%setup -q

%build
CFLAGS="-I%{_prefix}/include $RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix}"
if [ %{_lib} == 'lib64' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DLIB_SUFFIX=64"
fi
if [ %{zipext} == 'lzma' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DDOC_COMPRESS=lzma -DDOC_COMPRESS_EXT=lzma"
fi
if [ %{zipext} == 'xz' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DDOC_COMPRESS=xz -DDOC_COMPRESS_EXT=xz"
fi

cmake $CMAKE_FLAGS .
make
make info

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/bin/comap
%{_prefix}/bin/mica
%{_prefix}/share/info/comap.info.%{zipext}
%{_prefix}/share/man/man1/comap.1.%{zipext}
%{_prefix}/share/man/man1/mica.1.%{zipext}

%changelog
* Sat Apr 18 2015 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.5.2-1
- Extended output files
* Thu Oct 09 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.5.1-1
- New mapping procedure
- Simplified interface
- Bug fixes...
* Tue Feb 28 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.4.1-1
* Tue Mar 29 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.4.0-1

