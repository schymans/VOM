function initialseed(node, onCompleteCallback)
{
   var myobj = { label: "runmodel", id: "runmodel", href: "runmodel.html", target:"basefrm" };
   var tmpNode35 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "writepars", id: "writepars", href: "writepars.html", target:"basefrm" };
   var tmpNode66 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "writeloop", id: "writeloop", href: "writeloop.html", target:"basefrm" };
   var tmpNode67 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode67.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
