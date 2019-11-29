function modules__f90(node, onCompleteCallback)
{
   var myobj = { label: "vegwatmod", id: "vegwatmod", href: "vegwatmod.html", target:"basefrm" };
   var tmpNode4 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode4.isLeaf = true; 
   var myobj = { label: "vegmod", id: "vegmod", href: "vegmod.html", target:"basefrm" };
   var tmpNode5 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode5.isLeaf = true; 
   var myobj = { label: "watmod", id: "watmod", href: "watmod.html", target:"basefrm" };
   var tmpNode6 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode6.isLeaf = true; 
   var myobj = { label: "sce_mod", id: "sce_mod", href: "sce_mod.html", target:"basefrm" };
   var tmpNode7 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode7.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
