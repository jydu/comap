%define name CoMap
%define version 1.3.0
%define release 1
%define _prefix /usr/local

Summary: CoMap: a program for decting molecular co-evolution.
Name: %{name}
Version: %{version}
Release: %{release}
Vendor: Julien Dutheil
Source: http://download.gna.org/comap/%{name}-%{version}.tar.gz
License: CeCILL 2
Group: System Environment/Libraries
BuildRoot: %{_builddir}/%{name}-root
Packager: Julien Dutheil
Prefix: %{_prefix}
Requires: Bpp-Utils = 1.2.0
Requires: Bpp-NumCalc = 1.4.0
Requires: Bpp-Seq = 1.4.0
Requires: Bpp-Phyl = 1.5.0

%description
CoMap is a program for substitution mapping and molecular co-evolution analysis.

%prep
%setup -q

%build
CFLAGS="-I%{_prefix}/include $RPM_OPT_FLAGS" ./configure --prefix=%{_prefix}
make

%install
rm -rf $RPM_BUILD_ROOT
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS COPYING INSTALL NEWS README ChangeLog
%{_prefix}/bin/comap
%{_infodir}/comap.info

%changelog
* Wed Jan 23 2008 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- CoMap 1.3.0 release
* Tue Jul 19 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- CoMap 1.2.0 release
* Mon Jul 31 2006 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- CoMap 1.0.0 release
* Fri Nov 16 2005 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- First draft of the spec file

