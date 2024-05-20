;; The 'guix.scm' file for Pyretechnics, for use by 'guix shell'.
;; To isolate your development environment from the rest of your system, run:
;;
;; $ guix time-machine -C channels.scm -- shell --container --network --link-profile -S /usr/bin/env=bin/env --share=$HOME/.ssh

(use-modules
 ((gnu packages base)            #:select (coreutils which))
 ((gnu packages bash)            #:select (bash))
 ((gnu packages check)           #:select (python-pytest))
 ((gnu packages emacs)           #:select (emacs-minimal))
 ((gnu packages emacs-xyz)       #:select (emacs-htmlize))
 ((gnu packages geo)             #:select (gdal))
 ((gnu packages less)            #:select (less))
 ((gnu packages python-build)    #:select (python-hatchling))
 ((gnu packages python-xyz)      #:select (python-numpy python-rasterio))
 ((gnu packages ssh)             #:select (openssh))
 ((gnu packages version-control) #:select (git))
 ((guix build-system pyproject)  #:select (pyproject-build-system))
 ((guix gexp)                    #:select (local-file))
 ((guix git-download)            #:select (git-predicate))
 ((guix licenses)                #:select (epl2.0))
 ((guix packages)                #:select (package))
 ((guix utils)                   #:select (current-source-directory)))

(define vcs-file?
  ;; Return true if the given file is under version control.
  (or (git-predicate (current-source-directory))
      (const #t)))

(package
 (name "python-pyretechnics")
 (version "2024.5.7")
 (source (local-file "."
                     "pyretechnics-checkout"
                     #:recursive? #t
                     #:select?    vcs-file?))
 (build-system pyproject-build-system) ; includes python-toolchain
 ;; The 'build' phase runs these commands:
 ;;   import sys, importlib, json
 ;;   config_settings = json.loads (sys.argv[3])
 ;;   builder = importlib.import_module(sys.argv[1])
 ;;   builder.build_wheel(sys.argv[2], config_settings=config_settings)
 ;; The 'check' phase runs this command:
 ;;   pytest -vv
 ;; See file:/run/current-system/profile/share/guile/site/3.0/guix/build/pyproject-build-system.scm for more info.
 (arguments `(#:configure-flags '()
              #:test-backend    'pytest
              #:test-flags      '()))
 (native-inputs (list
                 ;; Shell utilities
                 bash
                 coreutils
                 which
                 less
                 ;; Version control
                 git
                 openssh
                 ;; GIS utilities
                 gdal
                 ;; Build tools
                 emacs-minimal
                 emacs-htmlize
                 python-hatchling
                 python-pytest))
 (propagated-inputs (list
                     ;; Python dependency libraries
                     python-numpy
                     python-rasterio))
 (synopsis "A Python library for simulating fire behavior in a variety of ways.")
 (description "A Python library for simulating fire behavior in a variety of ways.")
 (home-page "https://github.com/pyregence/pyretechnics/")
 (license epl2.0))
