%define _prefix /usr

URL: http://bioweb.me/comap/

Name: comap
Version: 1.5.5
Release: 1%{?dist}
License: CECILL-2.0
Vendor: Julien Dutheil
Source: %{name}-%{version}.tar.gz
Summary: The CoMap package
Group: Productivity/Scientific/Other

Requires: libbpp-phyl12 = 2.4.0
Requires: libbpp-seq12 = 2.4.0
Requires: libbpp-core' = 2.4.0

BuildRoot: %{_builddir}/%{name}-root
BuildRequires: cmake >= 2.8.11
BuildRequires: gcc-c++ >= 4.7.0
BuildRequires: groff
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core4 = 2.4.0
BuildRequires: libbpp-core-devel = 2.4.0
BuildRequires: libbpp-seq12 = 2.4.0
BuildRequires: libbpp-seq-devel = 2.4.0
BuildRequires: libbpp-phyl12 = 2.4.0
BuildRequires: libbpp-phyl-devel = 2.4.0


AutoReq: yes
AutoProv: yes
%if 0%{?mdkversion}
%if 0%{?mdkversion} >= 201100
BuildRequires: xz
%define compress_program xz
%else
BuildRequires: lzma
%define compress_program lzma
%endif
%else
BuildRequires: gzip
%define compress_program gzip
%endif

%description
Includes programs:
 - CoMap, (co)substitution mapping and coevolution detection,
 - MICA, Mutual Information Coevolution Analysis.

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DCOMPRESS_PROGRAM=%{compress_program}"
cmake $CMAKE_FLAGS .
make

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
%{_prefix}/share/info/comap.info.*
%{_prefix}/share/man/man1/comap.1.*
%{_prefix}/share/man/man1/mica.1.*

%changelog
* Thu Apr 19 2018 Julien Dutheil <dutheil@evolbio.mpg.de> 1.5.5-1
- Compatibility update with Bio++ 2.4.0.
- Migrated code to c++11.
- New mapping visualization tools (R).
* Thu Jun 08 2017 Julien Dutheil <dutheil@evolbio.mpg.de> 1.5.4-1
- Compatibility update with Bio++ 2.3.1.
* Sun May 21 2017 Julien Dutheil <dutheil@evolbio.mpg.de> 1.5.3-1
- New conditional group randomization procedure
- Compatibility update with Bio++ 2.3.0.
* Thu Oct 09 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.5.1-1
- New mapping procedure
- Simplified interface
- Bug fixes...
* Tue Feb 28 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.4.1-1
* Tue Mar 29 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.4.0-1

