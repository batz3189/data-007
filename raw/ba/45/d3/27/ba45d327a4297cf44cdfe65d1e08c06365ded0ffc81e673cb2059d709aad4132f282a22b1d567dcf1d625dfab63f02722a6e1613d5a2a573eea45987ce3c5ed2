pragma solidity ^0.4.0;
contract EnumDemo {
    enum ActionChoices{GoLeft,GoRight}
    
    ActionChoices public choice;
    
    ActionChoices defaultChoice = ActionChoices.GoRight;
    
    event ActionChoicesE (ActionChoices _ac,uint _u);
    function setChoice(){
        choice = ActionChoices.GoRight;
        //ActionChoicesE(choice,uint(choice));
    }
    
    function getChoice() returns(ActionChoices){
         return  choice;
          ActionChoicesE(choice,uint(choice));
    }
    
    function getDefault() returns(uint){
        return uint(defaultChoice);
    }
}