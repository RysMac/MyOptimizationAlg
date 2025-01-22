:Begin:
:Function: SMSSetLinkOption
:Pattern: SMSSetLinkOption["SubproblemAceGenTest",{i_Integer,j_Integer}]
:Arguments: {i,j}
:ArgumentTypes: {Integer,Integer}
:ReturnType:   Manual
:End:

:Begin:
:Function: SMSLinkNoEvaluations
:Pattern: SMSLinkNoEvaluations["SubproblemAceGenTest"]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType:   Manual
:End:

:Evaluate: SMSLinkNoEvaluations[]:=SMSLinkNoEvaluations[SMSSessionName];

:Evaluate: SMSSetLinkOptions[i__Rule]:=SMSSetLinkOptions[SMSSessionName,i];
:Evaluate: SMSSetLinkOptions[s_String,Rule[i_String,j_]]:=SMSSetLinkOption[s,{i,j}/.{{"SparseArray",True}->{0,2},{"SparseArray",False}->{0,1},{"SparseArray",Automatic}->{0,0},{"PauseOnExit",True}->{1,1},{"PauseOnExit",False}->{1,0},_:>(Print["Incorrect option: ",Rule[i,j]];Abort[])}];
:Evaluate: SMSSetLinkOptions[s_String,i__Rule]:=Map[SMSSetLinkOptions[s,#]&,{i}];

:Begin:
:Function: SubproblemAceGenTestMathLink
:Pattern: SubproblemAceGenTest[solution_?(ArrayQ[#,1,(Head[#]==Real || Head[#]==Integer&)] && Dimensions[#]==={2} &)]
:Arguments: {solution}
:ArgumentTypes: {Manual}
:ReturnType:   Manual
:End: