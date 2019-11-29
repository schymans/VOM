function vom_read_input(node, onCompleteCallback)
{
   var myobj = { label: "read_commandline", id: "read_commandline", href: "read_commandline.html", target:"basefrm" };
   var tmpNode7 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
