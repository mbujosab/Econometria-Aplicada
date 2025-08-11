%
% Author: Rolf Niepraschk, https://github.com/rolfn/kalenderRN
%
% inspired by http://www.texample.net/tikz/examples/a-calender-for-doublesided-din-a4/
% Original author: Robert Krause
% License: Creative Commons attribution license (CC-BY-SA)

%% Cambiar línea 333 (\sbox{\RN@pageII}{\RN@makeKalender{2019-01-01}{2019-06-30}}% )
%% donde el año es el segundo año del curso...

\setcounter{errorcontextlines}{100}
\NeedsTeXFormat{LaTeX2e}[1999/12/01]
\ProvidesClass{kalenderRN}[2016/07/25 v0.3b Calendar class (RN)]

\newif\ifRN@print \RN@printfalse

\DeclareOption{print}{\RN@printtrue}

\DeclareOption*{\PassOptionsToPackage{\CurrentOption}{xcolor}}
\ProcessOptions\relax
\LoadClass[a4paper,10pt]{article}
\ifRN@print
  \RequirePackage{graphicx}
  \PassOptionsToPackage{portrait}{geometry}
\else
  \PassOptionsToPackage{landscape}{geometry}
\fi
\RequirePackage[a4paper,margin={0mm,0mm}]{geometry}

\newif\ifRN@XeOrLua \RN@XeOrLuafalse
\RequirePackage{ifluatex,ifxetex}
\ifxetex\RN@XeOrLuatrue\else\ifluatex\RN@XeOrLuatrue\fi\fi

\ifRN@XeOrLua
  \RequirePackage{fontspec}
  \setsansfont[Ligatures=TeX]{TeX Gyre Heros}
\else
  \RequirePackage[T1]{fontenc}
  \RequirePackage{textcomp}
  \RequirePackage{lmodern}
  \RequirePackage{tgheros}
\fi

\renewcommand*\familydefault{\sfdefault}

\RequirePackage{etoolbox}
\AtEndPreamble{% "etoolbox"
  \ifRN@XeOrLua\else
    \@ifpackageloaded{inputenc}{}{%
      \RequirePackage[utf8]{inputenc}%
    }%
  \fi
}

%\RequirePackage[main=english]{babel}
\RequirePackage[main=spanish]{babel}

\RequirePackage[spanish]{translator}

\RequirePackage{tikz}
\usetikzlibrary{calc}
\usetikzlibrary{calendar}

\RequirePackage{varwidth}

\newcommand*\RN@year{}
\edef\RN@year{\the\numexpr\year}
\newcommand*\RN@yearS{}
\edef\RN@yearS{\the\numexpr\year+1}

% \makeatletter
% \let\c@nextyear\nextyear
% \makeatother
% \newcommand*{\thenextyear}{\Roman{year}}

% \newcommand*\RN@nextyear{}
% \edef\RN@nextyear{\the\numexpr\nextyear}

\newsavebox\RN@pageI
\newsavebox\RN@pageII
\newsavebox\RN@tempbox
\newdimen\RN@tempdima
\newdimen\RN@tempdimb
\newcommand*\RN@CalendarioEscolar@I{}
\newcommand*\RN@CalendarioEscolar@II{}

\newcommand*\RN@clases@I{}
\newcommand*\RN@clases@II{}

\newcommand*\RN@examens@I{}
\newcommand*\RN@examens@II{}

\newcommand*\RN@backgrounds{}
\newcommand*\RN@periods{}
\newcommand*\RN@title{}
\newcommand*\RN@tempa{}
\newcommand*\RN@tempb{}
\newcounter{weeknumber}
\newif\ifRN@showweeknumbers \RN@showweeknumbersfalse

\pgfkeys{%
  /CAL/.is family, /CAL,
  showweeknumbers/.is if = RN@showweeknumbers,
  title/.estore in = \RN@title,
  year/.estore in = \RN@year,
  default/.style={eventColor=medium,monthBGcolor=dark,monthColor=white,%
    workdayColor=none,saturdayColor=bright,sundayColor=medium,titleColor=dark}, %,examenColor=medium}, %,claseColor=medium
%  default/.style={claseColor=medium,monthBGcolor=lecture,monthColor=white,%
%    workdayColor=none,saturdayColor=bright,sundayColor=medium,titleColor=lecture},
%  default/.style={examenColor=medium,monthBGcolor=lecture,monthColor=white,%
%    workdayColor=none,saturdayColor=bright,sundayColor=medium,titleColor=examen},
  titleColor/.estore in = \RN@titleColor,
  eventColor/.estore in = \RN@eventColor,
  claseColor/.estore in = \RN@claseColor,
  examenColor/.estore in = \RN@examenColor,
  monthBGcolor/.estore in = \RN@monthBGcolor,
  monthColor/.estore in = \RN@monthColor,
  workdayColor/.estore in = \RN@workdayColor,
  saturdayColor/.estore in = \RN@saturdayColor,
  sundayColor/.estore in = \RN@sundayColor
}
\pgfkeys{%
  /CAL/.is family, /CAL,
  period/.cd,
  default/.style={color=period,name={}},
  color/.estore in = \RN@period@color,
  name/.estore in = \RN@period@name
}
\pgfkeys{%
  /CAL/.is family, /CAL,
  event/.cd,
  color/.estore in = \RN@eventColor
}
\pgfkeys{%
  /CAL/.is family, /CAL,
  clase/.cd,
  color/.estore in = \RN@claseColor
}
\pgfkeys{%
  /CAL/.is family, /CAL,
  examen/.cd,
  color/.estore in = \RN@examenColor
}

\newcommand*\period[2]{%
  \@ifnextchar[{\period@i{#1}{#2}}{\period@i{#1}{#2}[]}
}
\newcommand*\period@i{}
\def\period@i#1#2[#3]{%
  \pgfkeys{/CAL/period/.cd, default, #3}%
  \edef\@tempa{%
    \noexpand\ifdate{between=#1 and #2}{%
      \noexpand\tikzset{every day/.style={fill=\RN@period@color}}}{}%
  }%
  \expandafter\g@addto@macro\expandafter\RN@periods\expandafter{\@tempa}%
  \ifx\RN@period@name\@empty\else
    \pgfcalendar{}{#1}{#2}{%
      \edef\RN@tempa{\pgfcalendarcurrentyear-\pgfcalendarcurrentmonth-\pgfcalendarcurrentday}%
      \@ifundefined{event@\RN@tempa}{%
        \event{\RN@tempa}{\RN@period@name}%
        \edef\RN@tempb{\noexpand\expandafter\noexpand\gdef
          \noexpand\csname event@\RN@tempa\noexpand\endcsname{\RN@period@name}}%
        \RN@tempb
      }{}%
    }%
  \fi
}
\newcommand*\event{%
  \@ifstar{\@tempswatrue\RN@event@i}{\@tempswafalse\RN@event@i}%
}
\newcommand*\RN@event@i[2]{%
  \@ifnextchar[{\RN@event@ii{#1}{#2}}{\RN@event@ii{#1}{#2}[]}
}
\newcommand*\RN@event@ii{}%
\def\RN@event@ii#1#2[#3]{%
  \pgfkeys{/CAL/event/.cd, #3}%
  \@ifundefined{event@#1}{%
    \protected@edef\@tempa{\noexpand\@ifundefined{pgf@sh@ns@cal-#1}{}{%
      \noexpand\node[anchor=north west] at ($(cal-#1.north east)-(.6em,0em)$){%
        \noexpand\makebox[0pt][r]{\noexpand\resizebox{!}{1.15ex}{%
          \noexpand\varwidth{4.6cm}#2\noexpand\strut\noexpand\endvarwidth}}};}%
    }%
    \@namedef{event@#1}{#2}%
    \pgfcalendarifdate{#1}{at most=06-last}% first half-year?
    {\expandafter\g@addto@macro\expandafter\RN@CalendarioEscolar@I\expandafter{\@tempa}}%
    {\expandafter\g@addto@macro\expandafter\RN@CalendarioEscolar@II\expandafter{\@tempa}}%
  }{}%
  \if@tempswa
    \edef\@tempa{%
      \noexpand\ifdate{equals=#1}{\noexpand\tikzset{%
        every day/.style={fill=\RN@eventColor}}}{}
    }%
    \expandafter\g@addto@macro\expandafter\RN@backgrounds\expandafter{\@tempa}%
  \fi
}

\newcommand*\salaO{%
  \@ifstar{\@tempswatrue\RN@clase@i}{\@tempswafalse\RN@clase@i}%
}

\newcommand*\clase{%
  \stepcounter{Lecc}
  \@ifstar{\@tempswatrue\RN@clase@i}{\@tempswafalse\RN@clase@i}%
}

\newcommand*\claseYmedia{%
  \stepcounter{Lecc}
  \@ifstar{\@tempswatrue\RN@clase@i}{\@tempswafalse\RN@clase@i}%
}

\newcommand*\RN@clase@i[2]{%
  \@ifnextchar[{\RN@clase@ii{#1}{#2}}{\RN@clase@ii{#1}{#2}[]}
}
\newcommand*\RN@clase@ii{}%
\def\RN@clase@ii#1#2[#3]{%
  \pgfkeys{/CAL/clase/.cd, #3}%
  \@ifundefined{clase@#1}{%
    \protected@edef\@tempa{\noexpand\@ifundefined{pgf@sh@ns@cal-#1}{}{%
      \noexpand\node[anchor=north west] at ($(cal-#1.north east)-(.6em,0em)$){%
        \noexpand\makebox[0pt][r]{\noexpand\resizebox{!}{1.15ex}{%
          \noexpand\varwidth{4.6cm}#2\noexpand\strut\noexpand\endvarwidth}}};}%
    }%
    \@namedef{clase@#1}{#2}%
    \pgfcalendarifdate{#1}{at most=06-last}% first half-year?
    {\expandafter\g@addto@macro\expandafter\RN@clases@I\expandafter{\@tempa}}%
    {\expandafter\g@addto@macro\expandafter\RN@clases@II\expandafter{\@tempa}}%
  }{}%
  \if@tempswa
    \edef\@tempa{%
      \noexpand\ifdate{equals=#1}{\noexpand\tikzset{%
        every day/.style={fill=\RN@claseColor}}}{}
    }%
    \expandafter\g@addto@macro\expandafter\RN@backgrounds\expandafter{\@tempa}%
  \fi
}

\newcommand*\examen{%
  \@ifstar{\@tempswatrue\RN@examen@i}{\@tempswafalse\RN@examen@i}%
}
\newcommand*\RN@examen@i[2]{%
  \@ifnextchar[{\RN@examen@ii{#1}{#2}}{\RN@examen@ii{#1}{#2}[]}
}
\newcommand*\RN@examen@ii{}%
\def\RN@examen@ii#1#2[#3]{%
  \pgfkeys{/CAL/examen/.cd, #3}%
  \@ifundefined{examen@#1}{%
    \protected@edef\@tempa{\noexpand\@ifundefined{pgf@sh@ns@cal-#1}{}{%
      \noexpand\node[anchor=north west] at ($(cal-#1.north east)-(.6em,0em)$){%
        \noexpand\makebox[0pt][r]{\noexpand\resizebox{!}{1.15ex}{%
          \noexpand\varwidth{4.6cm}#2\noexpand\strut\noexpand\endvarwidth}}};}%
    }%
    \@namedef{examen@#1}{#2}%
    \pgfcalendarifdate{#1}{at most=06-last}% first half-year?
    {\expandafter\g@addto@macro\expandafter\RN@examens@I\expandafter{\@tempa}}%
    {\expandafter\g@addto@macro\expandafter\RN@examens@II\expandafter{\@tempa}}%
  }{}%
  \if@tempswa
    \edef\@tempa{%
      \noexpand\ifdate{equals=#1}{\noexpand\tikzset{%
        every day/.style={fill=\RN@examenColor}}}{}
    }%
    \expandafter\g@addto@macro\expandafter\RN@backgrounds\expandafter{\@tempa}%
  \fi
}

\definecolor{dark}{HTML}{A57C00}
\colorlet{medium}{dark!35}  % domingos
\colorlet{bright}{dark!20}  % sábados

\colorlet{period}{gray!30}  % no usado
\definecolor{examen}{HTML}{CC0000}
\colorlet{postex}{examen!20}
\colorlet{ex}{examen!50}
\colorlet{preex}{examen!20}
\definecolor{lecture}{HTML}{BF94E4}
\colorlet{lecc}{gray!35}
\definecolor{laboratorio}{HTML}{87A96B}
\definecolor{clase}{HTML}{848482} % días con lección

\newcommand\RN@@title{%
  \raisebox{-.3\height}{\parbox[b]{\wd\RN@tempbox}{%
    \fontsize{60}{70}\selectfont\color{\RN@titleColor}%
    %\textbf{\RN@year}~
    \textbf{\Huge\RN@title}\hfill
  }}%
}

\newcommand\RN@makeKalender[2]{%
  \sbox\RN@tempbox{%
    \begin{tikzpicture}[every day/.style={anchor=north}]%
      \calendar[dates=#1 to #2,name=cal,day text={\%d=},day yshift=3em,%]%
        day code={%
          \node[name=\pgfcalendarsuggestedname,every day,shape=rectangle,%
            minimum height=.53cm,text width=4.4cm,draw=black]{\tikzdaytext};%
            \draw(-1.8cm,-.1ex) node[anchor=west]{\footnotesize
              \pgfcalendarweekdayshortname{\pgfcalendarcurrentweekday}%
              \ifRN@showweeknumbers
                \ifdate{Monday}{~~\stepcounter{weeknumber}\normalsize[%
                  \mbox{\footnotesize\arabic{weeknumber}}]}{}%
              \fi
            };%
        },%
        execute before day scope={%
          \ifdate{day of month=1}{%
            \pgftransformxshift{4.8cm}%
            % Print month name
            \draw (0,0)node [shape=rectangle,minimum height=.53cm,%
              text width=4.4cm,fill=\RN@monthBGcolor,draw=\RN@monthBGcolor,%
                text=\RN@monthColor,text centered]{%
              \textbf{\pgfcalendarmonthname{\pgfcalendarcurrentmonth}~\pgfcalendarcurrentyear}};}{}%
          \ifdate{workday}{\tikzset{every day/.style={fill=\RN@workdayColor}}}{}%
          \RN@periods
          \ifdate{Saturday}{\tikzset{every day/.style={fill=\RN@saturdayColor}}}{}%
          \ifdate{Sunday}{\tikzset{every day/.style={fill=\RN@sundayColor}}}{}%
          \RN@backgrounds
        },%
        execute at begin day scope={%
          \pgftransformyshift{-.53*\pgfcalendarcurrentday cm}%
        }%
      ];%
      \pgfcalendarifdate{#1}{at most=06-last}{\RN@CalendarioEscolar@I}{\RN@CalendarioEscolar@II}%
      \pgfcalendarifdate{#1}{at most=06-last}{\RN@clases@I}{\RN@clases@II}%
      \pgfcalendarifdate{#1}{at most=06-last}{\RN@examens@I}{\RN@examens@II}%
    \end{tikzpicture}%
  }%
  \RN@tempdima=\ifRN@print\paperwidth\else\paperheight\fi
  \RN@tempdimb=\ifRN@print\paperheight\else\paperwidth\fi
  \parbox[c][\RN@tempdima]{\RN@tempdimb}{%
    \vfill
    \centering
    \RN@@title\par
    \vfill
    \usebox\RN@tempbox\par
    \vfill
  }%
}

\newcommand*\makeKalender[1][]{%
  \pgfkeys{/CAL/.cd, default, #1}%
  \setcounter{weeknumber}{0}%
  \begingroup
    \let\year=\RN@year
    \InputIfFileExists{clases.tex}{}{}%
    \InputIfFileExists{CalendarioEscolar.tex}{}{}%
  \endgroup
  \sbox{\RN@pageI}{\RN@makeKalender{\RN@year-09-01}{\RN@yearS-01-last}}%
  %\sbox{\RN@pageI}{\RN@makeKalender{\RN@year-09-01}{\RN@year-12-last}}%
  %\sbox{\RN@pageI}{\RN@makeKalender{\RN@year-07-01}{\RN@year-12-31}}%
  \sbox{\RN@pageII}{\RN@makeKalender{\RN@yearS-01-01}{\RN@yearS-02-last}}%       %% MBB
  \ifRN@print
    \rotatebox[origin=c]{90}{\usebox\RN@pageI}%
  \else
    \usebox\RN@pageI
  \fi
  %\newpage
  %\ifRN@print
  %  \rotatebox[origin=c]{-90}{\usebox\RN@pageII}%
  %\else
  %  \usebox\RN@pageII
  %\fi
}
\setlength\parindent{0pt}
\pagestyle{empty}
\endinput
