Get["evolucion_n1.m"];
Get["evolucion_n2.m"];
Get["evolucion_n3.m"];
Get["evolucion_n4.m"];
(*

In[40]:= NonlinearModelFit[evolucionn1, a x + c, {a,c},x][x]
Out[40]= 0.000179491 +1.63659 x
In[41]:= NonlinearModelFit[evolucionn2, a x + c, {a,c},x][x]
Out[41]= 0.000575063 +3.14502 x
In[48]:= 3.1450192661757406`/2
Out[48]= 1.57251
In[42]:= NonlinearModelFit[evolucionn3, a x + c, {a,c},x][x]
Out[42]= 0.00038385 +4.64635 x
In[49]:= 4.646345721446062`/3
Out[49]= 1.54878
In[43]:= NonlinearModelFit[evolucionn4, a x + c, {a,c},x][x]
Out[43]= -0.00265322+6.37662 x
In[50]:= 6.376621892764`/4
Out[50]= 1.59416


*)

Show[ListPlot[evolucionn1,AxesLabel->{"k","\[Lambda]"},PlotLabel->"Evolución del modo normal",PlotStyle->Thick], Plot[0.00017949118146448747` +1.6365913043527673` x,{x,0,0.01}]];

Show[ListPlot[evolucionn2,AxesLabel->{"k","\[Lambda]"},PlotLabel->"Evolución del modo normal",PlotStyle->Thick], Plot[0.0005750626354865649` +3.1450192661757406` x,{x,0,0.01}]]

Show[ListPlot[Delete[evolucionn3,1],AxesLabel->{"k","\[Lambda]"},PlotLabel->"Evolución del modo normal",PlotStyle->Thick], Plot[0.000383850309984158` +4.646345721446062` x,{x,0,0.01}]]

Show[ListPlot[evolucionn4, AxesLabel -> {"k", "\[Lambda]"}, 
  PlotLabel -> "Evolución del modo normal", PlotStyle -> Thick], 
 Plot[-0.002653217312251198` + 6.376621892764` x, {x, 0, 1/100}]]
