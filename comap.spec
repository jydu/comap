%define _basename comap
%define _version 1.4.1
%define _release 1
%define _prefix /usr

URL: http://home.gna.org/comap/

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: Julien Dutheil
Source: http://biopp.univ-montp2.fr/repos/sources/%{_basename}-%{_version}.tar.gz
Summary: The CoMap package
Group: Productivity/Scientific/Other

Requires: libbpp-phyl9 = 2.0.3
Requires: libbpp-seq9 = 2.0.3
Requires: libbpp-core2 = 2.0.3

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.6.0
BuildRequires: gcc-c++ >= 4.0.0
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core2 = 2.0.3
BuildRequires: libbpp-core-devel = 2.0.3
BuildRequires: libbpp-seq9 = 2.0.3
BuildRequires: libbpp-seq-devel = 2.0.3
BuildRequires: libbpp-phyl9 = 2.0.3
BuildRequires: libbpp-phyl-devel = 2.0.3

AutoReq: yes
AutoProv: yes
%if 0%{?mdkversion}
%define zipext xz
%else
%define zipext gz
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
* Tue Feb 28 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.4.1-1
* Tue Mar 29 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.4.0-1

