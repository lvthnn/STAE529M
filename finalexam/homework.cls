\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{homework}[2023/09/22 Homework LaTeX class]

\LoadClass{report}
\RequirePackage{fancyhdr}
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}
\RequirePackage{lastpage}
\RequirePackage{extramarks}
\RequirePackage[margin=1in]{geometry}
\RequirePackage{booktabs}
\RequirePackage{setspace}
\RequirePackage{enumitem}
\RequirePackage{parskip}
\RequirePackage{xcolor}
\RequirePackage{graphicx}
\RequirePackage{courier}
\RequirePackage{physics}
\RequirePackage{etoolbox}
\usepackage{float}
% \RequirePackage[capitalise]{cleveref}
\RequirePackage{lipsum}

\newcommand\crule[3][black]{\textcolor{#1}{\rule{#2}{#3}}}
\newcommand{\RomanNumeralCaps}[1]
    {\MakeUppercase{\romannumeral #1}}  

\RequirePackage{amsmath, amssymb} 
\RequirePackage{mathtools}
\RequirePackage{stmaryrd}
% \RequirePackage{minted}
\RequirePackage{arrayjobx}
\RequirePackage{tikz}
\usetikzlibrary{arrows.meta}

% LANGUAGE AND FONT CONFIGURATION OPTIONS
% Language options
\newarray\Language
\Language(1)={icelandic}
\Language(2)={english}

\newarray\hmwkDepartment
\hmwkDepartment(1)={Háskóli Íslands}
\hmwkDepartment(2)={University of Iceland}

\newarray\Pages
\Pages(1)={bls. \thepage\ af \pageref{LastPage}}
\Pages(2)={page \thepage\ of \pageref{LastPage}}

\newcommand{\configureLang}[1]{
    \checkLanguage(\hmwkLang)
    \ifdefstring{\cachedata}{icelandic}
    {\RequirePackage[icelandic]{datetime2}\RequirePackage[icelandic]{babel}}
    {\RequirePackage[en-GB]{datetime2}\RequirePackage[english]{babel}}
}

% Font setup

\newarray\Font
\Font(1)={default}
\Font(2)={times}
\Font(3)={palatino}
\Font(4)={fourier}

\newcommand{\temp}{}

\newcommand{\configureFont}[1]{
    \checkFont(\hmwkFont)
    \ifdefstring{\cachedata}{times}{\RequirePackage{times}}{}
}

\graphicspath{{img}}
\allowdisplaybreaks

%\input{setup/mathsetup}

% DOCUMENT SETUP
\topmargin=-0.45in
\evensidemargin=0in
\oddsidemargin=0in
\textwidth=6.5in
\textheight=9.0in
\headsep=0.25in
\headheight=15pt

% PAGE STYLE
\pagestyle{fancy}
\lhead{\hmwkDepartment(\hmwkLang)}
\chead{\hmwkClassNumber\ — \hmwkClass}
\rhead{\today}
\lfoot{\lastxmark}
\cfoot{}
\rfoot{\Pages(\hmwkLang)}
\renewcommand\headrulewidth{0.4pt}
\renewcommand\footrulewidth{0.4pt}

% TITLE PAGE SETUP
\title{\hmwkTitle \ - \hmwkClass}
\author{\hmwkAuthor \\ \texttt{\hmwkAuthEmail}}
\date{}

\thispagestyle{fancy}

