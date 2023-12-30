/* File: codegen.cc
 * ----------------
 * Implementation for the CodeGenerator class. The methods don't do anything
 * too fancy, mostly just create objects of the various Tac instruction
 * classes and append them to the list.
 */

#include "codegen.h"
#include <string.h>
#include "tac.h"
#include "mips.h"
#include "ast_decl.h"
#include "errors.h"
  
CodeGenerator::CodeGenerator()
{
  code = new List<Instruction*>();
  curGlobalOffset = 0;
}

char *CodeGenerator::NewLabel()
{
  static int nextLabelNum = 0;
  char temp[10];
  sprintf(temp, "_L%d", nextLabelNum++);
  return strdup(temp);
}


Location *CodeGenerator::GenTempVar()
{
  static int nextTempNum;
  char temp[10];
  Location *result = NULL;
  sprintf(temp, "_tmp%d", nextTempNum++);
  return GenLocalVariable(temp);
}

  
Location *CodeGenerator::GenLocalVariable(const char *varName)
{            
    curStackOffset -= VarSize;
    return new Location(fpRelative, curStackOffset+4,  varName);
}

Location *CodeGenerator::GenGlobalVariable(const char *varName)
{
    curGlobalOffset += VarSize;
    return new Location(gpRelative, curGlobalOffset -4, varName);
}


Location *CodeGenerator::GenLoadConstant(int value)
{
  Location *result = GenTempVar();
  code->Append(new LoadConstant(result, value));
  return result;
}

Location *CodeGenerator::GenLoadConstant(const char *s)
{
  Location *result = GenTempVar();
  code->Append(new LoadStringConstant(result, s));
  return result;
} 

Location *CodeGenerator::GenLoadLabel(const char *label)
{
  Location *result = GenTempVar();
  code->Append(new LoadLabel(result, label));
  return result;
} 


void CodeGenerator::GenAssign(Location *dst, Location *src)
{
  code->Append(new Assign(dst, src));
}


Location *CodeGenerator::GenLoad(Location *ref, int offset)
{
  Location *result = GenTempVar();
  code->Append(new Load(result, ref, offset));
  return result;
}

void CodeGenerator::GenStore(Location *dst,Location *src, int offset)
{
  code->Append(new Store(dst, src, offset));
}


Location *CodeGenerator::GenBinaryOp(const char *opName, Location *op1,
						     Location *op2)
{
  Location *result = GenTempVar();
  code->Append(new BinaryOp(BinaryOp::OpCodeForName(opName), result, op1, op2));
  return result;
}


void CodeGenerator::GenLabel(const char *label)
{
  code->Append(new Label(label));
}

void CodeGenerator::GenIfZ(Location *test, const char *label)
{
  code->Append(new IfZ(test, label));
}

void CodeGenerator::GenGoto(const char *label)
{
  code->Append(new Goto(label));
}

void CodeGenerator::GenReturn(Location *val)
{
  code->Append(new Return(val));
}


BeginFunc *CodeGenerator::GenBeginFunc(FnDecl *fn)
{
  BeginFunc *result = new BeginFunc;
  code->Append(insideFn = result);
  List<VarDecl*> *formals = fn->GetFormals();
  int start = OffsetToFirstParam;
  if (fn->IsMethodDecl()) start += VarSize;
  for (int i = 0; i < formals->NumElements(); i++)
    formals->Nth(i)->rtLoc = new Location(fpRelative, i*VarSize + start, formals->Nth(i)->GetName());
  curStackOffset = OffsetToFirstLocal;
  return result;
}

void CodeGenerator::GenEndFunc()
{
  code->Append(new EndFunc());
  insideFn->SetFrameSize(OffsetToFirstLocal-curStackOffset);
  insideFn = NULL;
}

void CodeGenerator::GenPushParam(Location *param)
{
  code->Append(new PushParam(param));
}

void CodeGenerator::GenPopParams(int numBytesOfParams)
{
  Assert(numBytesOfParams >= 0 && numBytesOfParams % VarSize == 0); // sanity check
  if (numBytesOfParams > 0)
    code->Append(new PopParams(numBytesOfParams));
}

Location *CodeGenerator::GenLCall(const char *label, bool fnHasReturnValue)
{
  Location *result = fnHasReturnValue ? GenTempVar() : NULL;
  code->Append(new LCall(label, result));
  return result;
}
  
Location *CodeGenerator::GenFunctionCall(const char *fnLabel, List<Location*> *args, bool hasReturnValue)
{
  for (int i = args->NumElements()-1; i >= 0; i--) // push params right to left
    GenPushParam(args->Nth(i));
  Location *result = GenLCall(fnLabel, hasReturnValue);
  GenPopParams(args->NumElements()*VarSize);
  return result;
}

Location *CodeGenerator::GenACall(Location *fnAddr, bool fnHasReturnValue)
{
  Location *result = fnHasReturnValue ? GenTempVar() : NULL;
  code->Append(new ACall(fnAddr, result));
  return result;
}
  
Location *CodeGenerator::GenMethodCall(Location *rcvr,
			     Location *meth, List<Location*> *args, bool fnHasReturnValue)
{
  for (int i = args->NumElements()-1; i >= 0; i--)
    GenPushParam(args->Nth(i));
  GenPushParam(rcvr);	// hidden "this" parameter
  Location *result= GenACall(meth, fnHasReturnValue);
  GenPopParams((args->NumElements()+1)*VarSize);
  return result;
}
 
 
static struct _builtin {
  const char *label;
  int numArgs;
  bool hasReturn;
} builtins[] =
 {{"_Alloc", 1, true},
  {"_ReadLine", 0, true},
  {"_ReadInteger", 0, true},
  {"_StringEqual", 2, true},
  {"_PrintInt", 1, false},
  {"_PrintString", 1, false},
  {"_PrintBool", 1, false},
  {"_Halt", 0, false}};

Location *CodeGenerator::GenBuiltInCall(BuiltIn bn,Location *arg1, Location *arg2)
{
  Assert(bn >= 0 && bn < NumBuiltIns);
  struct _builtin *b = &builtins[bn];
  Location *result = NULL;

  if (b->hasReturn) result = GenTempVar();
                // verify appropriate number of non-NULL arguments given
  Assert((b->numArgs == 0 && !arg1 && !arg2)
	|| (b->numArgs == 1 && arg1 && !arg2)
	|| (b->numArgs == 2 && arg1 && arg2));
  if (arg2) code->Append(new PushParam(arg2));
  if (arg1) code->Append(new PushParam(arg1));
  code->Append(new LCall(b->label, result));
  GenPopParams(VarSize*b->numArgs);
  return result;
}


void CodeGenerator::GenVTable(const char *className, List<const char *> *methodLabels)
{
  code->Append(new VTable(className, methodLabels));
}


void CodeGenerator::DoFinalCodeGen()
{
  createCFG();
  if (IsDebugOn("tac")) { // if debug don't translate to mips, just print Tac
    for (int i = 0; i < code->NumElements(); i++)
	code->Nth(i)->Print();
   }  else {
     Mips mips;
     mips.EmitPreamble();
     for (int i = 0; i < code->NumElements(); i++)
	 code->Nth(i)->Emit(&mips);
  }
}

void CodeGenerator::createCFG()
{
  Instruction *prev = NULL;
  //spot position of labels to jump to
  Hashtable<Instruction*> *labelTable = new Hashtable<Instruction*>;
  for (int i = 0; i < code->NumElements(); i++) {
    Label *label = dynamic_cast<Label*>(code->Nth(i));
    if (label) {
      labelTable->Enter(label->GetName(), code->Nth(i+1));
    }
   }

  //create CFG
  for (int i = 0; i < code->NumElements(); i++) {
    Instruction *instr = code->Nth(i);
    IfZ *ifz = dynamic_cast<IfZ*>(instr);
    Goto *gt = dynamic_cast<Goto*>(instr);
    if(dynamic_cast<EndFunc*>(instr)){
      continue;
    }else if(ifz){
      Instruction *next = labelTable->Lookup(ifz->GetLabel());
      ifz->successors.Append(next);
      ifz->successors.Append(code->Nth(i+1));
    }else if(gt){
      Instruction *next = labelTable->Lookup(gt->GetLabel());
      gt->successors.Append(next);
    }else{
      instr->successors.Append(code->Nth(i+1));
    }
  }

  liveVariableAnalysis();
}

void CodeGenerator::liveVariableAnalysis()
{
  bool changed = true;
  while(changed){
    changed = false;
    for(int i = code->NumElements()-1; i >= 0; i--){
      Instruction *instr = code->Nth(i);
      List<Location*> *out = new List<Location*>;
      for(int j = 0; j < instr->successors.NumElements(); j++){// OUT[TAC] = Union(IN[SUCC(TAC)]) 
        Instruction *succ = instr->successors.Nth(j);
        for(int k = 0; k < succ->liveIn->NumElements(); k++){
          out->Append(succ->liveIn->Nth(k));
        }
      }
      out->Unique();//remove duplicates
      if (out != instr->liveOut)
      {
        changed = true;
        instr->liveOut = out;
      }else{
        continue;
      }
      Location* k = instr->GetKill();
      List<Location*> *in = new List<Location*>; // IN'[TAC] = OUT[TAC] - KILL(TAC) + GEN(TAC) 
      for(int j = 0; j < out->NumElements(); j++){
        if((!k) || k!= out->Nth(j)){
          in->Append(out->Nth(j));
        }
      }
      for(int j = 0; j < instr->GetGen()->NumElements(); j++){
          in->Append(instr->GetGen()->Nth(j));
      }
      in->Unique();
      instr->liveIn = in;
      
    }
  }
  //print out live variables
  /*
  for(int i = 0; i < code->NumElements(); i++){
    Instruction *instr = code->Nth(i);
    printf("%s: ", instr->GetPrintNameForTac());
    for(int j = 0; j < instr->liveVars.NumElements(); j++){
      printf("%s ", instr->liveVars.Nth(j)->GetName());
    }
    printf("\n");
  }
  
  */
  
}

void CodeGenerator::BuildInterferenceGraph()
{
    InterferenceGraph_t* current = NULL;

    for (int i = 0; i < code->NumElements(); i++)
    {
        auto tac = code->Nth(i);
        if (auto beginfunc_tac = dynamic_cast<BeginFunc*> (tac))
            current = &(beginfunc_tac->interference_graph);

        if (current)
        {
            for (auto from_tac : *(tac->live_vars_in))
            {
                if (current->find(from_tac) == current->end())
                    (*current)[from_tac] = {};
                
                for (auto to_tac : *(tac->live_vars_in))
                {
                    if (from_tac != to_tac)
                        (*current)[from_tac].insert(to_tac);
                }
            }

            for (auto kill_tac : *(tac->GetKills()))
            {
                if (current->find(kill_tac) == current->end())
                    (*current)[kill_tac] = {};
                
                for (auto out_tac : *(tac->live_vars_out))
                {
                    if (kill_tac != out_tac)
                    {
                        (*current)[kill_tac].insert(out_tac);
                        (*current)[out_tac].insert(kill_tac);
                    }
                }
            }
        }
    }
}


Location *CodeGenerator::GenArrayLen(Location *array)
{
  return GenLoad(array, -4);
}

Location *CodeGenerator::GenNew(const char *vTableLabel, int instanceSize)
{
  Location *size = GenLoadConstant(instanceSize);
  Location *result = GenBuiltInCall(Alloc, size);
  Location *vt = GenLoadLabel(vTableLabel);
  GenStore(result, vt);
  return result;
}


Location *CodeGenerator::GenDynamicDispatch(Location *rcvr, int vtableOffset, List<Location*> *args, bool hasReturnValue)
{
  Location *vptr = GenLoad(rcvr); // load vptr
  Assert(vtableOffset >= 0);
  Location *m = GenLoad(vptr, vtableOffset*4);
  return GenMethodCall(rcvr, m, args, hasReturnValue);
}

// all variables (ints, bools, ptrs, arrays) are 4 bytes in for code generation
// so this simplifies the math for offsets
Location *CodeGenerator::GenSubscript(Location *array, Location *index)
{
  Location *zero = GenLoadConstant(0);
  Location *isNegative = GenBinaryOp("<", index, zero);
  Location *count = GenLoad(array, -4);
  Location *isWithinRange = GenBinaryOp("<", index, count);
  Location *pastEnd = GenBinaryOp("==", isWithinRange, zero);
  Location *outOfRange = GenBinaryOp("||", isNegative, pastEnd);
  const char *pastError = NewLabel();
  GenIfZ(outOfRange, pastError);
  GenHaltWithMessage(err_arr_out_of_bounds);
  GenLabel(pastError);
  Location *four = GenLoadConstant(VarSize);
  Location *offset = GenBinaryOp("*", four, index);
  Location *elem = GenBinaryOp("+", array, offset);
  return new Location(elem, 0); 
}



Location *CodeGenerator::GenNewArray(Location *numElems)
{
  Location *one = GenLoadConstant(1);
  Location *isNonpositive = GenBinaryOp("<", numElems, one);
  const char *pastError = NewLabel();
  GenIfZ(isNonpositive, pastError);
  GenHaltWithMessage(err_arr_bad_size);
  GenLabel(pastError);

  // need (numElems+1)*VarSize total bytes (extra 1 is for length)
  Location *arraySize = GenLoadConstant(1);
  Location *num = GenBinaryOp("+", arraySize, numElems);
  Location *four = GenLoadConstant(VarSize);
  Location *bytes = GenBinaryOp("*", num, four);
  Location *result = GenBuiltInCall(Alloc, bytes);
  GenStore(result, numElems);
  return GenBinaryOp("+", result, four);
}

void CodeGenerator::GenHaltWithMessage(const char *message)
{
   Location *msg = GenLoadConstant(message);
   GenBuiltInCall(PrintString, msg);
   GenBuiltInCall(Halt, NULL);
}
