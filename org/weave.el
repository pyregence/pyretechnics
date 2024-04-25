#!/usr/bin/env -S emacs -Q --script

;;==========================================================
;; Configure org's PDF export
;;==========================================================

(require 'org)

(setq org-latex-classes                '(("article" "\\documentclass[11pt]{article}"
                                          ("\\section{%s}"       . "\\section*{%s}")
                                          ("\\subsection{%s}"    . "\\subsection*{%s}")
                                          ("\\subsubsection{%s}" . "\\subsubsection*{%s}")
                                          ("\\paragraph{%s}"     . "\\paragraph*{%s}")
                                          ("\\subparagraph{%s}"  . "\\subparagraph*{%s}")))
      org-export-latex-hyperref-format "\\ref{%s}"
      org-latex-default-packages-alist '(("AUTO"     "inputenc"  t)
                                         (""         "graphicx"  t)
                                         (""         "longtable" nil)
                                         (""         "float"     nil)
                                         (""         "wrapfig"   nil)
                                         (""         "rotating"  nil)
                                         ("normalem" "ulem"      t)
                                         (""         "amsmath"   t)
                                         (""         "textcomp"  t)
                                         (""         "marvosym"  t)
                                         (""         "wasysym"   t)
                                         (""         "amssymb"   t)
                                         (""         "minted"    nil)
                                         (""         "hyperref"  nil)
                                         "\\tolerance=1000")
      org-latex-listings               'minted
      org-latex-minted-options         '(("linenos"         "false")
                                         ("frame"           "single")
                                         ("fontsize"        "\\scriptsize")
                                         ("baselinestretch" "1.0"))
      org-latex-pdf-process            '("pdflatex -shell-escape -interaction nonstopmode -output-directory %o %f"
                                         "bibtex %b"
                                         "pdflatex -shell-escape -interaction nonstopmode -output-directory %o %f"
                                         "pdflatex -shell-escape -interaction nonstopmode -output-directory %o %f")
      make-backup-files                nil)

;;==========================================================
;; Now load Pyretechnics.org and weave its prose into a PDF
;;==========================================================

(find-file "Pyretechnics.org")
(org-latex-export-to-pdf)

;;==========================================================
;; Obligatory calling shell protection
;;==========================================================

(setq argv nil)
