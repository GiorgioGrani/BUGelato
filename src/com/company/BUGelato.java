package com.company;

import javax.swing.*;

public class BUGelato {
    private JButton printParamsButton;
    private JButton Simulate;
    private JPanel field;
    private JButton exit;
    private JPasswordField passwordField1;
    private JTextArea output;

    public static void main(String[] args) {
        BUGelato bug = new BUGelato();
        JFrame frame = new JFrame("BUGelato");
        frame.setContentPane(bug.field);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);

        bug.exit.addActionListener((event) -> System.exit(0));

        String s="";
        bug.Simulate.addActionListener((event) -> BUGelato.runMain(bug) );



    }

    private static void runMain(BUGelato bug){
        //String s = Main.basic();
        //bug.output.setText(s);
    }
}
