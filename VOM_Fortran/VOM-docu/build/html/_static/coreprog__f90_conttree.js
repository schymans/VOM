function coreprog__f90(node, onCompleteCallback)
{
   var myobj = { label: "vom", id: "vom", href: "vom.html", target:"basefrm" };
   var tmpNode2 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode2.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
