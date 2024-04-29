;; The 'guix.scm' file for Pyretechnics, for use by 'guix shell'.
;; To isolate your development environment from the rest of your system, run:
;;
;; $ guix shell --container --network --link-profile -S /usr/bin/env=bin/env --share=$HOME/.ssh

(use-modules
 ((gnu packages base)            #:select (coreutils which))
 ((gnu packages emacs)           #:select (emacs-minimal))
 ((gnu packages emacs-xyz)       #:select (emacs-htmlize))
 ((gnu packages geo)             #:select (gdal))
 ((gnu packages less)            #:select (less))
 ((gnu packages python)          #:select (python))
 ((gnu packages python-xyz)      #:select (python-numpy python-rasterio))
 ((gnu packages ssh)             #:select (openssh))
 ((gnu packages version-control) #:select (git))
 ((guix build-system python)     #:select (python-build-system))
 ((guix licenses)                #:select (epl2.0))
 ((guix packages)                #:select (package)))

(package
 (name "python-pyretechnics")
 (version "2024.04.29")
 (source #f)
 (build-system python-build-system)
 (native-inputs (list
                 ;; Shell utilities
                 coreutils
                 which
                 less
                 ;; Version control
                 git
                 openssh
                 ;; GIS utilities
                 gdal
                 ;; Build tools
                 python
                 emacs-minimal
                 emacs-htmlize))
 (propagated-inputs (list
                     ;; Python dependency libraries
                     python-numpy
                     python-rasterio))
 (synopsis "A Python library for simulating fire behavior in a variety of ways.")
 (description "A Python library for simulating fire behavior in a variety of ways.")
 (home-page "https://github.com/pyregence/pyretechnics/")
 (license epl2.0))
