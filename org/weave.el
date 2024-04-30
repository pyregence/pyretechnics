#!/usr/bin/env -S emacs -Q --script

;;==========================================================
;; Configure org's HTML export
;;==========================================================

(require 'org)
(require 'ox-html)
(require 'htmlize)
(require 'python)

(setq python-indent-guess-indent-offset nil)

;;==========================================================
;; Now load pyretechnics.org and weave it into an HTML file
;;==========================================================

(message "Exporting literate documentation to ../doc/pyretechnics.html")
(find-file "pyretechnics.org")
(org-html-export-to-html)

;;==========================================================
;; Obligatory calling shell protection
;;==========================================================

(setq argv nil)
