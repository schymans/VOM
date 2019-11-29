function sce_init(node, onCompleteCallback)
{
   var myobj = { label: "read_shufflevar", id: "read_shufflevar", href: "read_shufflevar.html", target:"basefrm" };
   var tmpNode17 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "open_output", id: "open_output", href: "open_output.html", target:"basefrm" };
   var tmpNode18 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
