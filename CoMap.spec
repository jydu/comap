%define name Coe
%define version 0.2.1
%define release 1
%define _prefix /usr/local

Summary: Coe: a program for decting molecular co-evolution.
Name: %{name}
Version: %{version}
Release: %{release}
Vendor: The Bio++ Project
Source: http://kimura.univ-montp2.fr/BioPP/Download/Sources/%{name}-%{version}.tar.gz
License: CeCILL 2
Group: System Environment/Libraries
BuildRoot: %{_builddir}/%{name}-root
Packager: Julien Dutheil
Prefix: %{_prefix}
Requires: Bpp-Utils = 1.0
Requires: Bpp-NumCalc = 1.0
Requires: Bpp-Seq = 1.0
Requires: Bpp-Phyl = 1.0
Requires: Coevolution = 0.1

%description
This is an unofficial program for molecular co-evolution analysis.
It is still under development and is provided for convenience.

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
%{_prefix}/bin/coe

%changelog
* Fri Nov 16 2005 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- First draft of the spec file

