function read_commandline(node, onCompleteCallback)
{
   var myobj = { label: "command_argument_count", id: "command_argument_count", href: "", target:"basefrm" };
   var tmpNode3 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode3.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
