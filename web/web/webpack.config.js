var path = require('path');
var webpack = require('webpack');

var isProd = (process.env.NODE_ENV === 'production');

module.exports = {
  entry: './src/main.js',

  output: {
    path: 'static/',
    filename: 'bundle.js'
  },

  module: {
    loaders: [
      {
        test: /\.jsx?$/,
        loader: 'babel-loader',
        include: [
          path.resolve(__dirname, "src")
        ],
        query: {
          presets: ['es2015', 'react']
        }
      }
    ]
  },

  plugins: isProd ? [
    new webpack.DefinePlugin({
      "process.env": { NODE_ENV: '"production"' }
    }),
    new webpack.optimize.UglifyJsPlugin({
      compress: { warnings: false }
    })
  ] : []
};
