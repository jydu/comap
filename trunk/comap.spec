%define name CoMap
%define version 1.4.1
%define release 1
%define _prefix /usr

Summary: The CoMap package.
Name: %{name}
Version: %{version}
Release: %{release}
Vendor: Julien Dutheil
Source: http://download.gna.org/comap/%{name}-%{version}.tar.gz
License: CeCILL 2
Group: Productivity/Scientific/Other
BuildRoot: %{_builddir}/%{name}-root
Packager: Julien Dutheil
Prefix: %{_prefix}
AutoReq: yes
AutoProv: yes

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
rm -rf $RPM_BUILD_ROOT
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
%{_prefix}/share/info/comap.info.gz
%{_prefix}/share/man/man1/comap.1.gz
%{_prefix}/share/man/man1/mica.1.gz

%changelog
* Tue Feb 28 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- CoMap 1.4.1 release
* Tue Mar 29 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- CoMap 1.4.0 release

