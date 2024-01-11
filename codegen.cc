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
#include <stack>
  
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
  BeginFunc *result = new BeginFunc();
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
  Mips* mips = interferenceGraph();
  if (IsDebugOn("tac")) { // if debug don't translate to mips, just print Tac
    for (int i = 0; i < code->NumElements(); i++)
	code->Nth(i)->Print();
   }  else {
     mips->EmitPreamble();
     for (int i = 0; i < code->NumElements(); i++)
	      code->Nth(i)->Emit(mips);
  }
}

void CodeGenerator::createCFG()
{
  //std::cout<<"createCFG"<<std::endl;
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
      if(i< code->NumElements()-1) ifz->successors.Append(code->Nth(i+1));
    }else if(gt){
      Instruction *next = labelTable->Lookup(gt->GetLabel());
      gt->successors.Append(next);
    }else{
      if(i< code->NumElements()-1) instr->successors.Append(code->Nth(i+1));
    }
  }
  //std::cout<<"liveVariableAnalysis"<<std::endl;
  liveVariableAnalysis();
}

bool CodeGenerator::compareLocList(List<Location*> *l1, List<Location*> *l2){
  if(l1->NumElements() != l2->NumElements()){
    return false;
  }
  for(int i = 0; i < l1->NumElements(); i++){
    //std::cout<<"compareLocList: "<<l1->Nth(i)->GetName()<<" "<<l2->Nth(i)->GetName()<<std::endl;
    if(l1->Nth(i)->GetName() != l2->Nth(i)->GetName() || l1->Nth(i)->GetSegment() != l2->Nth(i)->GetSegment() || l1->Nth(i)->GetOffset() != l2->Nth(i)->GetOffset()){
      return false;
    }
  }
  return true;
}

bool CodeGenerator::compareLoc(Location *l1, Location *l2){
  if( l1->GetName() == l2->GetName() && l1->GetSegment() == l2->GetSegment() && l1->GetOffset() == l2->GetOffset()){
    return true;
  }
  return false;
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
      instr->liveOut = out;
      Location* k = instr->GetKill();
      List<Location*>* g = instr->GetGen();
      List<Location*> *in = new List<Location*>; // IN'[TAC] = OUT[TAC] - KILL(TAC) + GEN(TAC) 
      //std :: cout << "in: " << instr->liveIn->NumElements() << std :: endl;
      for(int j = 0; j < out->NumElements(); j++){
        if((!k) || out->Nth(j)!=k){
        
          in->Append(out->Nth(j));
        }
      }
      for(int j = 0; j < g->NumElements(); j++){
          in->Append(g->Nth(j));
      }
      in->Unique();
      if (!in && !instr->liveIn)
      {
        continue;
      }else if ((!in && instr->liveIn) || (in && !instr->liveIn))
      {
        changed = true;
        instr->liveIn = in;
      }else if (!compareLocList(instr->liveIn, in))
      { 
        changed = true;
        instr->liveIn = in;
      }else{
        continue;
      }
      
      
    }
  }
  //print out live variables
   /*debugging
  for(int i = 0; i < code->NumElements(); i++){
    Instruction *instr = code->Nth(i);
    for(int j = 0; j < instr->liveOut->NumElements(); j++){
      printf("%s ", instr->liveOut->Nth(j)->GetName());
    }
    printf("\n");
  }
 
  */
  

  
}

List<Location*>*  CodeGenerator::reorderLoc(List<Location*> *l){
  List<Location*> *newList = new List<Location*>;
  //std::stack<Location*> *stack = new std::stack<Location*>;
  int max = 0;
  Location *maxLoc = nullptr;
  while (l->NumElements() > 0)
  {
    for(int i = 0; i < l->NumElements(); i++){
      Location *loc = l->Nth(i);
      if (loc->GetInterferences()->NumElements() > max)
      {
        max = loc->GetInterferences()->NumElements();
        maxLoc = loc;
      }
    }

    newList->Append(maxLoc);
    l->Remove(maxLoc);
    max = 0;
  }
  
  return newList;
  
}

Mips* CodeGenerator::interferenceGraph(){
  List<Location*>* graphlst = nullptr;
  Mips* mips = new Mips;
  //std::stack<Hashtable<Hashtable<Location*>*>*> graphStack;
  for (size_t i = 0; i < code->NumElements(); i++)
  {
    Instruction *instr = code->Nth(i);
    Label *label = dynamic_cast<Label*>(instr);
    BeginFunc *begin = NULL;
    if(i<code->NumElements()-1){
      begin = dynamic_cast<BeginFunc*>(code->Nth(i+1));
    }

    if (label && begin)
    {
      //std :: cout << "begin" << std :: endl;
      graphlst = begin->graphList;
    }
    if (instr->liveOut->NumElements() == 0){
      if (dynamic_cast<EndFunc*>(instr))
      {
        graphlst=reorderLoc(graphlst);
        mips->AllocateRegisters(graphlst);
      }
      continue;
    }
    //std::cout<<"interferenceGraph: "<<instr->liveOut->NumElements()<<std::endl;
    //kill
    Location *kill = instr->GetKill();
    if (kill)
    {
      
      for (size_t j = 0; j < instr->liveOut->NumElements(); j++)
      {
        //std::cout<<"interferenceGraph: "<<instr->liveOut->Nth(j)->GetName()<<std::endl;
        Location *loc = instr->liveOut->Nth(j);
        kill->AddInterference(loc);
        loc->AddInterference(kill);
      }
        graphlst->Remove(kill);
        graphlst->Append(kill);
      
      
    }

    //other
    Location *loc ;
    for (size_t j = 0; j < instr->liveOut->NumElements(); j++)
    {
      loc = instr->liveOut->Nth(j);
      //std :: cout << "interferenceGraph: " << loc->GetName() << std :: endl;
      for (size_t k = 0; k < instr->liveOut->NumElements(); k++)
      {
        if (j == k)
        {
          continue;
        }
        //std :: cout << "addi: " << instr->liveOut->Nth(k)->GetName() << std :: endl;
        loc->AddInterference(instr->liveOut->Nth(k));
      }
      
      if(kill) loc->AddInterference(kill);
      
      graphlst->Remove(loc);
      graphlst->Append(loc);
      //std :: cout << "loc had int: " << loc->GetName() << loc->GetInterferences()->NumElements() <<graphlst->NumElements()<< std :: endl;
    }//also enter oneself, avoid empty hashtable when removing
  }

  //print out interference graph
 /*
  for (size_t i = 0; i < code->NumElements(); i++)
  {
    Instruction *instr = code->Nth(i);
    for (size_t j = 0; j < instr->liveOut->NumElements(); j++)
    {
      Location *loc = instr->liveOut->Nth(j);
      printf("%s: ", loc->GetName());
      for (size_t k = 0; k < loc->GetInterferences()->NumElements(); k++)
      {
        printf("%s ", loc->GetInterferences()->Nth(k)->GetName());
      }
      printf("\n");
    }
  } 
*/
  return mips;
  
  
}

/*
void CodeGenerator::kColoring(Hashtable<Hashtable<Location*>*>* curGraph){
  Iterator<Hashtable<Location*>*> iter = curGraph->GetIterator();
  Hashtable<Location*> *kv;
  List<Location*> *nodeList = new List<Location*>;
  //int maxDegree = 0;
  //Hashtable<Location*> *maxDegreeNode = nullptr;
  while ((kv = iter.GetNextValue()) != nullptr)//find max degree
  {
      Iterator<Location*> iter2 = kv->GetIterator();
      Location *kv2;
      while ((kv2 = iter2.GetNextValue()) != nullptr)
      {
        
      }
      
      
  }

}

*/

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
