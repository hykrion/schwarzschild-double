#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" ${1+"$@"}

package require Tk

# -----------------------------------------------------------------------------
# @brief  Interfaz gráfica para el caso de Schwarzschild
#
#         Simplemente escribimos los parámetros en modo texto en el fichero
#         parameters.txt
# -----------------------------------------------------------------------------

# Datos habituales para los cálculos
set uiCommon(RangoInf) -60.0
set uiCommon(RangoSup) 60.0
set uiCommon(l) 0.0
set uiCommon(rS) 2.0
set uiCommon(wMin) 0.02
set uiCommon(wMax) 1.0
set uiCommon(nW) 100
set uiCommon(nL) 5
# Datos habituales para la gráfica de la onda
set uiCommon(X_LOW) -60.0
set uiCommon(X_HIG) 60.0
set uiCommon(Y_LOW) -2.0
set uiCommon(Y_HIG) 2.0

# -----------------------------------------------------------------------------
# @brief  Escribir configuración en ficheros gnuplot
# -----------------------------------------------------------------------------
proc writeGnuplotFile {fileName fileContent} {
  set fd [open $fileName w+]
  
  foreach line $fileContent {
    puts $fd $line
  }
  
  close $fd
}

# -----------------------------------------------------------------------------
# @brief  Modificar algunos parámetros de gnuplot. Reescribimos los x.gnu
# -----------------------------------------------------------------------------
proc writeGnuplotSettings {} {
  # Onda real / imaginaria
  lappend fileContent {set title 'Schwarzschild'}
  lappend fileContent "set xrange \[$::uiCommon(X_LOW) : $::uiCommon(X_HIG)\]"
  lappend fileContent "set yrange \[$::uiCommon(Y_LOW) : $::uiCommon(Y_HIG)\]"
  lappend fileContent {plot 'wave.txt' using 1:2 with lines lc "green" title "Φ Re", 'wave.txt' using 1:3 with lines lc "red" title "Φ Im", 'wave.txt' using 1:4 with lines title "V(x)"}
  
  writeGnuplotFile "plot.gnu" $fileContent

  # Tortuga
  set fileContent {}
  lappend fileContent {set title 'Coordenadas tortuga'}
  lappend fileContent {set xrange [-60.0 : 60.0]}
  lappend fileContent {set yrange [-0.1 : 60.0]}
  lappend fileContent {plot 'turtle.txt' using 1:2 with lines lc "blue" title "x[y] calculado", 'turtle.txt' using 3:1 with lines lc "red" title "y[x] ana.inv.", 'turtle.txt' using 4:1 with lines lc "green" title "y[x] cal.inv."}

  writeGnuplotFile "plot-turtle.gnu" $fileContent

  # Potencial
  set fileContent {}
  lappend fileContent {set title 'Potencial'}
  lappend fileContent {set xrange [-60.0 : 60.0]}
  lappend fileContent {set yrange [-0.05 : 1.2]}
  lappend fileContent {plot 'wave.txt' using 1:4 with lines lc "blue" title "V(x)"}
  
  writeGnuplotFile "plot-v.gnu" $fileContent
  
  # Coeficientes RT
  set fileContent {}
  lappend fileContent {set title 'Coeficientes R y T'}
  lappend fileContent {set xrange [0.0 : 1.0]}
  lappend fileContent {set yrange [0.0 : 1.2]}
  lappend fileContent {plot 'coefficients.txt' using 1:2 with lines lc "green" title "R", 'coefficients.txt' using 1:3 with lines lc "red" title "T", 'coefficients.txt' using 1:4 with lines title "R + T"}
  
  writeGnuplotFile "plot-k.gnu" $fileContent
  
  # Sigma
  set fileContent {}
  lappend fileContent {set title 'σl'}
  lappend fileContent {set xrange [0.0 : 1.0]}
  lappend fileContent {set yrange [0.0 : 120.0]}
  lappend fileContent {plot 'sigma-l.txt' using 1:2 with lines lc "red" title "σl / w"}
  
  writeGnuplotFile "plot-sigma.gnu" $fileContent
}

# -----------------------------------------------------------------------------
# @brief  Ejecutables con opciones
# -----------------------------------------------------------------------------
proc runExe {fileNameAndOption {delay 3000}} {
  after $delay
  exec {*}$fileNameAndOption &
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Comandos
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# -----------------------------------------------------------------------------
# @brief  Salvar los parámetros el el fichero 'parameters.txt', realizar los
#         cálculos mediante el ejecutable y ejecutar gnuplot para visualizar los
#         datos.
#         El número y orden han de coincidir con los declarados en 'ui.c'.
# -----------------------------------------------------------------------------
proc cmd_visualizar {} {
  # Salvar los parámetros
  set fd [open parameters.txt w+]
  
  puts $fd $::uiCommon(RangoInf); #0
  puts $fd $::uiCommon(RangoSup); #1
  puts $fd $::uiCommon(l);        #2
  puts $fd $::uiCommon(rS);       #3
  puts $fd $::uiCommon(wMin);     #4
  puts $fd $::uiCommon(wMax);     #5
  puts $fd $::uiCommon(nW);       #6
  puts $fd $::uiCommon(nL);       #7
  
  close $fd
  
  # Realizar los cálculos
  set command [list schwarzschild.exe]
  exec {*}$command
  after 1000

  writeGnuplotSettings
  cmd_replot
}

# -----------------------------------------------------------------------------
# @brief  Redibujar las gráficas
# -----------------------------------------------------------------------------
proc cmd_replot {} {
  # Coordenada tortuga
  set gnuplotCmd [list gnuplot.exe -p plot-turtle.gnu]
  runExe $gnuplotCmd

  # Dibujar los resultados del potencial y la función de onda
  set gnuplotCmd [list gnuplot.exe -p plot.gnu]
  runExe $gnuplotCmd
  
  set gnuplotCmd [list gnuplot.exe -p plot-v.gnu]
  runExe $gnuplotCmd

  # Coeficientes RT
  set gnuplotCmd [list gnuplot.exe -p plot-k.gnu]
  runExe $gnuplotCmd

  # Sigma
  set gnuplotCmd [list gnuplot.exe -p plot-sigma.gnu]
  runExe $gnuplotCmd
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# UI
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# -----------------------------------------------------------------------------
# @brief  Aquí vamos a modificar los parámetros del método numérico
# -----------------------------------------------------------------------------
proc ui_numerico {} {
  set mainFrame [frame .ntb.frmNumerico.frmMain]

  set w1 [label $mainFrame.lblRango -text "Rango: "]
  set w2 [label $mainFrame.lblRangoInf -text "Inf"]
  set w3 [entry $mainFrame.entRangoInf -textvariable uiCommon(RangoInf) -width 5 -background white]
  set w4 [label $mainFrame.lblRangoSup -text "Sup"]
  set w5 [entry $mainFrame.entRangoSup -textvariable uiCommon(RangoSup) -width 5 -background white]
  grid $w1 $w2 $w3 $w4 $w5
    
  set w1 [button $mainFrame.btnVisualizar -text "Visualizar" -command cmd_visualizar]
  set w2 [button $mainFrame.btnRePlot -text "Re Plot" -command cmd_replot]
  grid $w1 $w2
  
  # Layout
  grid $mainFrame -sticky ew
}

# -----------------------------------------------------------------------------
# @brief  Modificar algunos aspectos de la gráfica mostrada por gnuplot
# -----------------------------------------------------------------------------
proc ui_fisico {} {
  set mainFrame [frame .ntb.frmFisico.frmMain]
  
  foreach id {l nL rS wMin wMax nW} {
    set w1 [label $mainFrame.lbl$id -text "$id: "]
    set w2 [entry $mainFrame.ent$id -textvariable uiCommon($id) -width 5 -background white]
    
    grid $w1 $w2
  }
  
  # Layout
  grid $mainFrame -sticky ew
}

# -----------------------------------------------------------------------------
# @brief  Aquí vamos a modificar los parámetros físicos
# -----------------------------------------------------------------------------
proc ui_gnuplot {} {
  set mainFrame [frame .ntb.frmGnuplot.frmMain]
  
  foreach id {X_LOW X_HIG Y_LOW Y_HIG} {
    set w1 [label $mainFrame.lbl$id -text "$id: "]
    set w2 [entry $mainFrame.ent$id -textvariable uiCommon($id) -width 5 -background white]
    
    grid $w1 $w2
  }
  
  # Layout
  grid $mainFrame -sticky ew
}

# -----------------------------------------------------------------------------
# @brief  Iniciar el sencillo GUI
# -----------------------------------------------------------------------------
proc ui_init {} {
  ttk::notebook .ntb
  .ntb add [frame .ntb.frmNumerico] -text "Num\u00E9rico"
  .ntb add [frame .ntb.frmFisico] -text "F\u00EDsico"
  .ntb add [frame .ntb.frmGnuplot] -text "gnuplot"
  .ntb select .ntb.frmNumerico
  
  ui_numerico
  ui_fisico
  ui_gnuplot

  grid .ntb
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Main
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc main {} {
  ui_init
}

main
