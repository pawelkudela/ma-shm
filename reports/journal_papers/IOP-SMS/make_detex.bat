SETLOCAL
set original="Parallel-SMS-R1.tex"
::set output="l:\Praca\Grammarly-LaTeX-man\grammarly_web\detex-out.txt"
set output="\\tsclient\E\grammarly_web\detex-out.txt"
::set detex="l:\Praca\Grammarly-LaTeX-man\opendetex\detex.exe"
set detex=%~dp0\..\..\..\bin\external\opendetex\detex.exe
%detex% -l %original%>%output%
ENDLOCAL
