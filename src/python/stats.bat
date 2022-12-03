set R="C:\Program Files\R\R-4.1.2\bin\R.exe"

REM stats
%R% --vanilla < fxtn.R > fxtn.out
%R% --vanilla < wm.R > wm.out
%R% --vanilla --args 6 < tm.R > tm.out

